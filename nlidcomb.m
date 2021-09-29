function [params_est_v, x0_est_v, diag_data] = nlidcomb(one_sample, fs, u_data_V, y_data_V, ...
    y_states_v, y_mapping, param_lb_v, param_guess_v, param_ub_v, ...
    lbx_x0_free_v, x0_v, ubx_x0_free_v, x0_free_v, scale_v, N_group_size, ...
    lbx_X_guess_H, X_guess_H, ubx_X_guess_H, additional_constraints_lbg_v, additional_constraints_g, ...
    additional_constraints_ubg_v, enable_gaps_v, fopts)
% Partially constrained multiple shooting (PCMS) and prediction error method (PEM) for grey-box system identification.
% For N_group_size = 1: equivalent to multiple shooting (MS).
% For N_group_size = size(y_data_V,1): equivalent to single shooting (SS).
% For PEM use fit_against_X_guess_H_v and fix_states_to_X_guess_H_v.

import casadi.*
% _v: vertical vector 
% _h: horizontal vector
% _V: tall and skinny matrix (1st dimension longer)
% _H: long and fat matrix (2nd dimension longer)
% (some of these variables can be a character vector as well)
% X_: variables storing how the states evolve in time (in symbolic or numeric form)
% y_: variables storing how the outputs evolve in time (in symbolic or numeric form)
% u_: variables storing how the inputs evolve in time (in numeric form)
% x0_: the initial state from which we start the simulation
% solver_x0_v: this is a different x0, it is the initial values for all the decision variables for the solver


%% Handle fopts
if(~isfield(fopts,'fit_against_X_guess_H_v')), fopts.fit_against_X_guess_H_v = []; end
if(~isfield(fopts,'fix_states_to_X_guess_H_v')), fopts.fix_states_to_X_guess_H_v = []; end
if(~isfield(fopts,'disable_gaps_per_state_v')), fopts.disable_gaps_per_state_v = []; end
if(~isfield(fopts,'gaps_abs_max')), fopts.gaps_abs_max = 0; end
if(~isfield(fopts,'opti_mode')), fopts.opti_mode = false; end
if(~isfield(fopts,'ipopt')), fopts.ipopt = struct; end

%% Handle free/fixed x0
assert(all(size(x0_free_v) == size(x0_v)));
assert(size(x0_free_v,2) == 1);
x0_decvars_v = MX.sym('x0_decvars_v',sum(x0_free_v),1);
x0_free_indexes_v = find(x0_free_v==true);
x0_sym_v = MX(x0_v); %this is just numbers so far (in MX form), not symbols... Need to convert from double so that 
    %we can do the following assignment to individual elements in the matrix:
x0_sym_v(x0_free_indexes_v) = x0_decvars_v;
if isempty(lbx_x0_free_v), lbx_x0_free_v = -inf(size(x0_decvars_v)); end
if isempty(ubx_x0_free_v), ubx_x0_free_v = inf(size(x0_decvars_v)); end
assert(all(size(lbx_x0_free_v)==size(x0_decvars_v)))
assert(all(size(ubx_x0_free_v)==size(x0_decvars_v)))

%% Initialize parameters
%N_group_size: the number of steps we take within the group. The number of items to LS against can either be N_group_size or N_group_size+1 (for first group).
N_samples = size(u_data_V,1); %The number of input/output samples that we are fitting against.
assert(N_group_size > 0)
assert(N_samples > 1)
N_sample_groups = floor((N_samples-1)/N_group_size); %The number of groups that we'll .map() over
N_grouped_samples = N_group_size * N_sample_groups; %The number of samples in the groups that we'll LS against.  
if N_sample_groups > 0, N_grouped_samples = N_grouped_samples + 1; end %...so this shall be increased by 1 for x0.
N_remainder_samples = N_samples - N_grouped_samples; %The samples that are not in a group and we'll LS against them.
N_remainder_steps = N_samples - N_group_size * N_sample_groups - 1; %The number of steps we have to take in the remainder part.
disp(['nlidcomb: N_samples = ' num2str(N_samples) ', N_group_size = ' num2str(N_group_size) newline ...
    '   => N_sample_groups = ' num2str(N_sample_groups) ', N_remainder_samples = ' num2str(N_remainder_samples) ',' newline ...
    '      N_grouped_samples = ' num2str(N_grouped_samples) ', N_remainder_steps = ' num2str(N_remainder_steps)]);
% The reason for (N_samples-1) here is that we have an extra sample from x0 to fit against.
% We have to take that into consideration when calculating these values. The following outputs are expected:
% N_group_size = 3, N_samples = 11 -> N_sample_groups = 3, N_remainder_samples = 1
% N_group_size = 3, N_samples = 10 -> N_sample_groups = 3, N_remainder_samples = 0
% N_group_size = 3, N_samples = 9  -> N_sample_groups = 2, N_remainder_samples = 2
% N_group_size = 3, N_samples = 8  -> N_sample_groups = 2, N_remainder_samples = 1
% N_group_size = 3, N_samples = 6  -> N_sample_groups = 1, N_remainder_samples = 2
% N_group_size = 3, N_samples = 4  -> N_sample_groups = 1, N_remainder_samples = 0
% If there is no group, then all samples should go into the remainder:
% N_group_size = 3, N_samples = 3  -> N_sample_groups = 0, N_remainder_samples = 3 
% N_group_size = 3, N_samples = 2  -> N_sample_groups = 0, N_remainder_samples = 2  

N_states = size(x0_v,1);
N_free_states = N_states - size(fopts.fix_states_to_X_guess_H_v, 1);
assert(N_free_states >= 0);
dont_fix_states_to_X_guess_H_v = setdiff(1:N_states, fopts.fix_states_to_X_guess_H_v).'; %the states that we don't want to be fixed to X_guess_H

%N_threads = 2^floor(log2(feature('numcores'))); % If we want to measure timing very accurately, 
%we should always have 2^x threads to eliminate artifacts in graphs, thus we always round the number of cores down 
%to the nearest power of 2, e.g. from 6 we round down to 4.
N_threads = feature('numcores'); % This is for performance, not for a timing benchmark.
params_v = MX.sym('params_v',length(param_guess_v),1);
assert(size(param_guess_v,2) == 1)
assert(all(size(param_lb_v) == size(param_guess_v)))
assert(all(size(param_ub_v) == size(param_guess_v)))
assert(all(size(scale_v)    == size(param_guess_v)))
assert(all(size(u_data_V,1) == size(y_data_V,1)))

%% X_guess_H is handled in this part
if strcmp(X_guess_H, 'guess_mechatronic')
    yd_V = diff(y_data_V(1:end,y_states_v))*fs; 
    X_guess_H = [y_data_V(2:end,y_states_v) yd_V].'; %X_guess_H shall actually be an extension of x0_sym_v... But it should not contain it neither.
    if(size(X_guess_H,1)<N_states) 
        X_guess_H = [X_guess_H;zeros(N_states-size(X_guess_H,1),size(X_guess_H,2))];
    end
elseif strcmp(X_guess_H, 'guess_exact')
    X_guess_H = [y_data_V(2:end,y_states_v)].'; %TODO check this case
else %in this case we were given an initial approximation of how all the states do evolve, our only task is to filter and keep the relevant samples:
    X_guess_H = X_guess_H(:,1:(size(y_data_V,1)-1)); %if too many samples were provided by the user, we cut them off first.
end

%bounds empty? Let's make a matrix out of -inf the same size as X_guess_H
if isempty(lbx_X_guess_H), lbx_X_guess_H = -inf(size(X_guess_H)); end 
%bounds empty? Let's make a matrix out of inf the same size as X_guess_H
if isempty(ubx_X_guess_H), ubx_X_guess_H =  inf(size(X_guess_H)); end 
%bounds given as a single number? Let's make this the same size as X_guess_H
if all(size(lbx_X_guess_H)==[1 1]), lbx_X_guess_H = lbx_X_guess_H.*ones(size(X_guess_H)); end
if all(size(ubx_X_guess_H)==[1 1]), ubx_X_guess_H = ubx_X_guess_H.*ones(size(X_guess_H)); end
assert(all(size(lbx_X_guess_H)==size(X_guess_H)))
assert(all(size(ubx_X_guess_H)==size(X_guess_H)))

% Let's keep only one guess per gap
X_guess_full_H = X_guess_H;
X_guess_H = X_guess_H(:,N_group_size:N_group_size:end-1);
X_guess_free_states_H = X_guess_H(dont_fix_states_to_X_guess_H_v, :);
lbx_X_guess_H = lbx_X_guess_H(dont_fix_states_to_X_guess_H_v,N_group_size:N_group_size:end-1);
ubx_X_guess_H = ubx_X_guess_H(dont_fix_states_to_X_guess_H_v,N_group_size:N_group_size:end-1);

if strcmp(class(X_guess_free_states_H),'double'), X_guess_free_states_H = casadi.DM(X_guess_free_states_H); end
if strcmp(class(lbx_X_guess_H),'double'), lbx_X_guess_H = casadi.DM(lbx_X_guess_H); end
if strcmp(class(ubx_X_guess_H),'double'), ubx_X_guess_H = casadi.DM(ubx_X_guess_H); end

%% Create maps and mapaccums
mapaccum_opts = struct;
%mapaccum_opts.base = 8;
if N_sample_groups>0
    if N_group_size == 1
        % If we only have to take 1 step, then we use one_step directly instead of through mapaccum: 
        % this is more effective for MXFunction (but the same for SXFunction).
        sample_group_mapaccum = one_sample;
    else %if N_group_size > 1 (it cannot be <= 0 by assert)
        sample_group_mapaccum = one_sample.mapaccum('sample_group_mapaccum', N_group_size, mapaccum_opts); %The steps to be taken in a group always equal to N_group_size.
    end
    if N_sample_groups == 1
        sample_group_map_intermediate = sample_group_mapaccum;
    else %if N_sample_groups > 1 (we already know that it's > 0)
        sample_group_map_intermediate = sample_group_mapaccum.map(N_sample_groups, 'thread', N_threads);
    end
    sample_group_map = sample_group_map_intermediate;
    %if fopts.fix_states_to_X_guess_H_v
    %    sgm_states_H = sample_group_map_intermediate.mx_in{1};
    %    sgm_controls_H = sample_group_map_intermediate.mx_in{2};
    %    sgm_params_H = sample_group_map_intermediate.mx_in{3};
    %    sgm_states_partfixed_map_H = sgm_states_H; %so these are the states that are partially fixed
    %    sgm_states_partfixed_map_H(fopts.fix_states_to_X_guess_H_v,:) = X_guess_H(fopts.fix_states_to_X_guess_H_v,1:N_sample_groups);
    %    sgm_states_notfixed_H = sgm_states_H(dont_fix_states_to_X_guess_H_v,:);
    %    sample_group_map = Function('sample_group_map', ...
    %        {sgm_states_notfixed_H, sgm_controls_H, sgm_params_H}, ...
    %        {sample_group_map_intermediate(sgm_states_partfixed_map_H, sgm_controls_H, sgm_params_H)})
    %else
    %    sample_group_map = sample_group_map_intermediate;
    %end
end

if N_remainder_steps > 1
    remainder_samples_mapaccum_intermediate = one_sample.mapaccum('remainder_samples_mapaccum_intermediate', N_remainder_steps, mapaccum_opts);
elseif N_remainder_steps == 1 
    remainder_samples_mapaccum_intermediate = one_sample;
else %if N_remainder_steps == 0
    remainder_samples_mapaccum_intermediate = [];
end
remainder_samples_mapaccum = remainder_samples_mapaccum_intermediate;

%if remainder_samples_mapaccum && fopts.fix_states_to_X_guess_H_v && 
%    rsma_states_H = remainder_samples_mapaccum_intermediate.mx_in{1};
%    rsma_controls_H = remainder_samples_mapaccum_intermediate.mx_in{2};
%    rsma_params_H = remainder_samples_mapaccum_intermediate.mx_in{3};
%    rsma_states_partfixed_map_H = rsma_states_H; %so these are the states that are partially fixed
%    assert(size(X_guess_H,2)==N_sample_groups+1) %we check if X_guess_H has otherwise been used up for the groups correctly 
%    rsma_states_partfixed_map_H(fopts.fix_states_to_X_guess_H_v,:) = X_guess_H(fopts.fix_states_to_X_guess_H_v,end); %only the last item is kept normally
%    rsma_states_notfixed_H = rsma_states_H(dont_fix_states_to_X_guess_H_v,:);
%    remainder_samples_mapaccum = Function('remainder_samples_mapaccum', ...
%        {rsma_states_notfixed_H, rsma_controls_H, rsma_params_H}, ...
%        {remainder_samples_mapaccum_intermediate(rsma_states_partfixed_map_H, rsma_controls_H, rsma_params_H)})
%else
%    remainder_samples_mapaccum = remainder_samples_mapaccum_intermediate;
%end


if N_sample_groups
    X_grouped_simstart_free_states_without_x0_H = MX.sym('X_grouped_simstart_free_states_without_x0_H', N_free_states, N_sample_groups-1);
    X_grouped_simstart_without_x0_H = MX.zeros(N_states, N_sample_groups-1); %does not have x0_sym_v
    if(~isempty(fopts.fix_states_to_X_guess_H_v))
        X_grouped_simstart_without_x0_H(fopts.fix_states_to_X_guess_H_v,:) = X_guess_H(fopts.fix_states_to_X_guess_H_v,1:N_sample_groups-1);
    end
    if(~isempty(dont_fix_states_to_X_guess_H_v))
        X_grouped_simstart_without_x0_H(dont_fix_states_to_X_guess_H_v,:) = X_grouped_simstart_free_states_without_x0_H;
    end
    X_grouped_simstart_H = [x0_sym_v X_grouped_simstart_without_x0_H];
    X_grouped_simend_without_x0_H = sample_group_map(X_grouped_simstart_H, u_data_V(1:N_grouped_samples-1,:).', repmat(params_v.*scale_v,1,N_grouped_samples-1)); %TODO add asserts to all map and mapaccum so that the inputs sizes match the inputs of the function
else
    X_grouped_simstart_free_states_without_x0_H = MX.zeros(N_free_states, 0);
    X_grouped_simstart_without_x0_H = MX.zeros(N_states, 0);
    X_grouped_simstart_H = MX.zeros(N_states, 0);
    X_grouped_simend_without_x0_H = MX.zeros(N_states, 0);
end

if N_remainder_samples > 0 
    if N_sample_groups == 0 %handle the case when there are no groups --> it is equivalent to single shooting
        X_remainder_simstart_free_states_without_x0_v = MX.zeros(N_free_states, 0);
        X_remainder_simstart_without_x0_v = MX.zeros(N_states, 0);
        X_remainder_simstart_v = x0_sym_v;
    else
        X_remainder_simstart_free_states_without_x0_v = MX.sym('X_remainder_simstart_without_x0_v', N_free_states, 1);
        X_remainder_simstart_without_x0_v = MX.zeros(N_states, 1);
        assert(size(X_guess_H,2)==N_sample_groups) %we check if X_guess_H has otherwise been used up for the groups correctly. X_guess_H should have 1 column per sample group, - 1 column for x0 (it does not have) + 1 column for remainder = N_sample_groups columns
        if(~isempty(fopts.fix_states_to_X_guess_H_v))
            X_remainder_simstart_without_x0_v(fopts.fix_states_to_X_guess_H_v,:) = X_guess_H(fopts.fix_states_to_X_guess_H_v,end); %only the last item is kept normally
        end
        if(~isempty(dont_fix_states_to_X_guess_H_v))
            X_remainder_simstart_without_x0_v(dont_fix_states_to_X_guess_H_v,:) = X_remainder_simstart_free_states_without_x0_v;
        end
        X_remainder_simstart_v = X_remainder_simstart_without_x0_v;
    end
    if N_sample_groups == 0 %handle the case when there are no groups
        remainder_starting_index = 1;
    else
        remainder_starting_index = N_grouped_samples;
    end
    u_data_remainder_H = u_data_V(remainder_starting_index:end-1,:).';
    assert(size(u_data_remainder_H,2)==N_remainder_steps);
    X_remainder_simend_without_x0_H = remainder_samples_mapaccum(X_remainder_simstart_v.', u_data_remainder_H, repmat(params_v.*scale_v,1,N_remainder_steps)); %we cannot use the last u_data 
    % The thing X_remainder_simend_without_x0_H is without x0_sym_v means that it does not contain it as the first state as is.
    % However, it still contains the simulation one sample after it... So it actually depends on x0_sym_v.
else
    X_remainder_simstart_v = MX.zeros(N_states, 0);
    X_remainder_simend_without_x0_H = MX.zeros(N_states, 0);
    X_remainder_simstart_without_x0_v = MX.zeros(N_states, 0);
    X_remainder_simstart_free_states_without_x0_v = MX.zeros(N_free_states, 0);
end

X_simstart_without_x0_H = [X_grouped_simstart_without_x0_H X_remainder_simstart_without_x0_v];
X_simstart_free_states_without_x0_H = [X_grouped_simstart_free_states_without_x0_H X_remainder_simstart_free_states_without_x0_v];
X_simend_without_x0_H = [X_grouped_simend_without_x0_H X_remainder_simend_without_x0_H]; 
X_sim_H = [x0_sym_v X_simend_without_x0_H];
assert(size(X_sim_H,2)==size(u_data_V,1))
gap_states_v = sort(setdiff(dont_fix_states_to_X_guess_H_v,fopts.disable_gaps_per_state_v));
gaps_H = X_sim_H(gap_states_v,N_group_size+1:N_group_size:end-1)-X_simstart_without_x0_H(gap_states_v,:);
if ~(ischar(enable_gaps_v) && strcmp(enable_gaps_v,'all')) 
    if(isempty(enable_gaps_v)) % enable_gaps_v = [] will disable all gaps
        gaps_H = MX(size(gap_states_v,1), 0);
    else
        gaps_H = gaps_H(:,enable_gaps_v);
    end
    
end

if isempty(y_mapping)
    y_sim_V = X_sim_H(y_states_v,:).';
else
    y_mapping_map = y_mapping.map(N_samples,'thread',N_threads);
    y_sim_V = y_mapping_map(X_sim_H,u_data_V.',repmat(params_v.*scale_v,1,N_samples)).';
end

assert(all(size(X_guess_H)==size(X_simstart_without_x0_H)));
assert(all(size(X_guess_free_states_H)==size(X_simstart_free_states_without_x0_H)));

%the error that we want to minimize --> will go into the objective
if isempty(fopts.fit_against_X_guess_H_v)
    e_V = y_data_V-y_sim_V;
else
    e_V = X_sim_H(fopts.fit_against_X_guess_H_v,:) - [x0_v(fopts.fit_against_X_guess_H_v) X_guess_full_H(fopts.fit_against_X_guess_H_v,:)];
end

%the decision variables:
V_v = veccat(params_v, x0_decvars_v, X_simstart_free_states_without_x0_H(:));

%assemble constraints: gaps + additional constraints
g_v = gaps_H(:);
lbg_v = -ones(size(g_v))*fopts.gaps_abs_max;
ubg_v = ones(size(g_v))*fopts.gaps_abs_max;
if ~isempty(additional_constraints_g) 
    acg_H = additional_constraints_g(X_sim_H, u_data_V.', repmat(params_v.*scale_v,1,N_samples));
    g_v=veccat(g_v, acg_H(:)); %this corresponds to the whole state
    lbg_v = [lbg_v;additional_constraints_lbg_v(:)];
    ubg_v = [ubg_v;additional_constraints_ubg_v(:)];
end

%the objective:
obj = 0.5*sumsqr(e_V); 
%obj = det(1/size(e_V,1).*(e_V.'*e_V));
nlp = struct('x',V_v,'f',obj);
if size(g_v,1)~=0, nlp.g=g_v; end %add constraints if there are gaps/additional_constraints

solver_x0_v = veccat(param_guess_v, x0_v(x0_free_indexes_v), X_guess_free_states_H(:));
solver_x0_lbx = veccat(param_lb_v, lbx_x0_free_v, lbx_X_guess_H(:));
solver_x0_ubx = veccat(param_ub_v, ubx_x0_free_v, ubx_X_guess_H(:));

%For debugging:
%f_x0_evalf = evalf(substitute(nlp.f,nlp.x,veccat(solver_x0_v)))
%g_x0_evalf = evalf(substitute(nlp.g,nlp.x,veccat(solver_x0_v)))
%f_diag_data_evalf = evalf(substitute(nlp.f,nlp.x,veccat(casadi.DM([0.14842;0.1313;0.2798;0.2798]),diag_data.X(:,1:N_group_size:end))))
%g_diag_data_evalf = evalf(substitute(nlp.g,nlp.x,veccat(casadi.DM([0.14842;0.1313;0.2798;0.2798]),diag_data.X(:,1:N_group_size:end))))

if fopts.opti_mode
    opti = Opti();
    theta_opti_v = opti.variable(size(nlp.x,1),1);
    mk_opti_f = Function('mk_opti_f', {nlp.x},{nlp.f},{'x'},{'f'});
    opti.minimize(mk_opti_f(theta_opti_v));
    opti.subject_to(solver_x0_lbx <= theta_opti_v <= solver_x0_ubx);
    if(isfield(nlp,'g'))
        mk_opti_g = Function('mk_opti_g', {nlp.x},{nlp.g},{'x'},{'g'}); 
        opti.subject_to(lbg_v <= mk_opti_g(theta_opti_v) <= ubg_v)
        %opti.subject_to(mk_opti_g(theta_opti_v) == 0)
    end
    opti.solver('ipopt',struct,fopts.ipopt);
    opti.set_initial(theta_opti_v, solver_x0_v);
    opti_sol = opti.solve();
    sol.x = opti_sol.value(theta_opti_v); %So there we are, we have emulated nlpsol and the result `sol` also emulates it, but we have opti.debug! Yeah.
else
    solver = sysid_gauss_newton(e_V,nlp,V_v,fopts);
    sol = solver('x0',solver_x0_v.','lbg',lbg_v.','ubg',ubg_v.','lbx',solver_x0_lbx.','ubx',solver_x0_ubx.'); 
end
params_est_v = full(sol.x(1:size(param_guess_v,1))).*scale_v;
x0_est_v = x0_v;
x0_est_v(x0_free_indexes_v) = full(sol.x(size(param_guess_v,1)+(1:size(x0_decvars_v,1))));

diag_data = struct;
if ~fopts.opti_mode
    diag_data.solver = solver;
end
diag_data.size_X_guess_H = size(X_guess_H);
diag_data.size_X_guess_free_states_H = size(X_guess_free_states_H);
diag_data.size_X_grouped_simstart_H = size(X_grouped_simstart_H);
diag_data.size_X_grouped_simend_without_x0_H = size(X_grouped_simend_without_x0_H);
diag_data.size_X_remainder_simstart_without_x0_v = size(X_remainder_simstart_without_x0_v);
diag_data.size_X_remainder_simend_without_x0_H = size(X_remainder_simend_without_x0_H);
diag_data.X_simstart_free_states_without_x0_H = reshape(full(sol.x((size(param_guess_v,1)+size(x0_decvars_v,1)+1):end)),size(X_guess_free_states_H));
diag_data.X_guess_H = full(X_guess_H);
diag_data.X_guess_full_H = X_guess_full_H;
diag_data.X_simstart_without_x0_H = full(evalf(substitute(X_simstart_without_x0_H, nlp.x, sol.x)));
diag_data.size_solver_x0_v = size(solver_x0_v);
diag_data.size_gaps_H = size(gaps_H);
diag_data.size_y_sim_V = size(y_sim_V);
diag_data.X_sim_H = full(evalf(substitute(X_sim_H, nlp.x, sol.x)));
diag_data.N_free_states = N_free_states;
diag_data.N_sample_groups = N_sample_groups;
diag_data.N_grouped_samples = N_grouped_samples;
diag_data.N_remainder_samples = N_remainder_samples;