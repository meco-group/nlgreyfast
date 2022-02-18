%% Sweep study for the AMC2022 paper
% Running nlgreyfast on: 
% - 100 cases of different input signals, 
% - different methods including SS, MS, PCMS, PEM, PEM+SS
% Before running this, you need to MEX the file `emps_plant_atan_m` using MATLAB coder.
% Test call for MATLAB Coder: emps_plant_atan_m(0, [1;2], 2, 3, 4, 5, 6, {}) 
% For variable types in MATLAB Coder, you need to specify a 0x0 cell inside an 1x1 cell for varargin.

%% Initialize model parameters 
nlid_emps_init
stable_initpoint = false; %this option allows us to set whether we want the identification to start from a stable initialization point or an unstable initialization point

%% Generate input signals
% We are generating many different multisines of the same parameters but random phases -> u_data_v
% We are simulating the model with all of them -> y_data_v
N_experiment_y = 100;
u_data_v = {};
X_sim_H = {};
for data_i = 1:(N_experiment_y+1) % the additional data will be for validation
    N_block = 5000;
    N_block_repeat = 2;
    N = N_block * N_block_repeat;
    
    mfms_opts = mkfullms_opts(N_block,fs,1,4);
    [~,multisine_cf,u_data_v{data_i}] = mkfullms(mkfullms_random(mfms_opts),[],mfms_opts);
    u_data_v{data_i} = 900.*repmat(u_data_v{data_i}, N_block_repeat, 1);
    
    x0_v{data_i} = min(max(0.3*randn(2,1),-3),3);
    [~,X_sim_H{data_i}] = nlid_emps_sim_plant(x0_v{data_i}, u_data_v{data_i}, param_truth);
    y_data_v{data_i} = full(X_sim_H{data_i}(1,:)).';
end

%% Identification for all cases and methods
delete nlid_emps_sweep_logs.txt; diary nlid_emps_sweep_logs.txt; disp('Processing started:'); disp(datetime('now')); diary off

N_hybrid_cases = ceil(log2(N));
N_additional_methods = 3; %in addition to PCMS, we will add 3 additional cases: PEM, PEM+SS, nlgreyest
hybrid_group_sizes=2.^((N_hybrid_cases+N_additional_methods):-1:0); 
y_states = [1];

for thread_index = 1:(N_experiment_y*length(hybrid_group_sizes)) 
    initial_params_index=floor((thread_index-1)/length(hybrid_group_sizes))+1 %it's just the data number, not initial_params_index anymore
    hybrid_group_size_index=mod((thread_index-1),length(hybrid_group_sizes))+1
    hybrid_group_size = hybrid_group_sizes(hybrid_group_size_index);
    nlid_emps_ubx = inf(4,1);
    nlid_emps_lbx = -inf(4,1);
    fopts = struct;
    fopts.compiler_flags = '';
    fopts.ipopt.max_iter = 1000;
    x0_free_v = [false;false];
    param_guess_v = iif(stable_initpoint,[ 100;100;100;-100 ],[ 0.2;-10;0.2;10 ]);
    whatisit="D/nlgreyest/PEM/PEM+SS, unstable init, no time limits, multi excitation, no noise";
    tic;
    if hybrid_group_size_index == 1 %PEM
        logs{thread_index}.case_description = 'PEM';
        fopts = struct;
        fopts.fit_against_X_guess_H_v = [1; 2];
        fopts.fix_states_to_X_guess_H_v = [1; 2];
        pem_group_size = 1;
        [param_est_comb, x0_est_v, diag_comb] = nlidcomb(plant_next, fs, ...
            u_data_v{initial_params_index}, y_data_v{initial_params_index}, ...
            y_states, [], nlid_emps_lbx, param_guess_v, nlid_emps_ubx, ...
            [], x0_v{initial_params_index}, [], x0_free_v, scale, pem_group_size, ...
            [], 'guess_mechatronic', [],  [], [], ...
            [], 'all', fopts);
        logs{thread_index}.dcstats = diag_comb.solver.stats;
    elseif hybrid_group_size_index == 2 %PEM+SS
        logs{thread_index}.case_description = 'PEM2';
        %PEM
        fopts = struct;
        fopts.fit_against_X_guess_H_v = [1; 2];
        fopts.fix_states_to_X_guess_H_v = [1; 2];
        pem_group_size = 1;
        [param_est_st1, x0_est_st1_v, diag_st1] = nlidcomb(plant_next, fs, ...
            u_data_v{initial_params_index}, y_data_v{initial_params_index}, ...
            y_states, [], nlid_emps_lbx, param_guess_v, nlid_emps_ubx, ...
            [], x0_v{initial_params_index}, [], x0_free_v, scale, pem_group_size, ...
            [], 'guess_mechatronic', [],  [], [], ...
            [], 'all', fopts);
        %SS
        ss_group_size = size(u_data_v{initial_params_index},1);
        fopts = struct;
        [param_est_comb, x0_est_v, diag_st2] = nlidcomb(plant_next, fs, ...
            u_data_v{initial_params_index}, y_data_v{initial_params_index}, ...
            y_states, [], nlid_emps_lbx, param_est_st1, nlid_emps_ubx, ...
            [], x0_est_st1_v, [], x0_free_v, scale, ss_group_size, ...
            [], 'guess_mechatronic', [],  [], [], ...
            [], 'all', fopts);
        logs{thread_index}.dcstats_pem1 = diag_st1.solver.stats;
        logs{thread_index}.dcstats = diag_st2.solver.stats;
    elseif hybrid_group_size_index == 3 %nlgreyest
        logs{thread_index}.case_description = 'nlgreyest';
        fopts = struct;
        fopts.x0_free_v = x0_free_v;
        fopts.u_data_v = u_data_v{initial_params_index};
        fopts.y_data_v = y_data_v{initial_params_index};
        logs{thread_index}.dcstats = struct;
        try
            [param_est_comb, x0_est_v, nlgr_nlgreyest] = nlid_emps_nlgreyest('emps_plant_atan_m_mex', param_guess_v, x0_v{initial_params_index}, fopts);
            logs{thread_index}.nlgr_nlgreyest = nlgr_nlgreyest;
            logs{thread_index}.dcstats.return_status = 'Solver_Succeeded';
        catch ME
            if strcmp(ME.identifier,'Ident:estimation:infeasibleSimulation')||strcmp(ME.identifier,'Ident:general:NoiseVarianceNotReal')
                param_est_comb = param_guess_v;
                x0_est_v = x0_v{initial_params_index};
                logs{thread_index}.dcstats.return_status = ['Nlgreyest_Fail:' ME.identifier];
                logs{thread_index}.ME = ME;
            else
                rethrow(ME)
            end
        end
    elseif hybrid_group_size_index > N_additional_methods
        logs{thread_index}.case_description = ['D = ' num2str(hybrid_group_size)];
        [param_est_comb, x0_est_v, diag_comb] = nlidcomb(plant_next, fs, ...
            u_data_v{initial_params_index}, y_data_v{initial_params_index}, ...
            y_states, [], nlid_emps_lbx, param_guess_v, nlid_emps_ubx, ...
            [], x0_v{initial_params_index}, [], x0_free_v, scale, hybrid_group_size, ...
            [], 'guess_mechatronic', [],  [], [], ...
            [], 'all', fopts);
        dcstats = diag_comb.solver.stats;
        logs{thread_index}.time_casadi = dcstats.t_wall_nlp_f + dcstats.t_wall_nlp_g + dcstats.t_wall_nlp_grad + dcstats.t_wall_nlp_grad_f + dcstats.t_wall_nlp_hess_l + dcstats.t_wall_nlp_jac_g;
        logs{thread_index}.time_ipopt = dcstats.t_wall_total - logs{thread_index}.time_casadi;
        logs{thread_index}.dcstats = dcstats;
    else
        error('hybrid_group_size_index is invalid')
    end
    logs{thread_index}.time_toc = toc;
    rmserror = rms(param_est_comb - param_truth);
    logs{thread_index}.whatisit = whatisit;
    logs{thread_index}.initial_params_index = initial_params_index;
    logs{thread_index}.hybrid_group_size = hybrid_group_size;
    logs{thread_index}.rmserror = rmserror;
    logs{thread_index}.param_est_comb = param_est_comb;
    logs{thread_index}.lbx = nlid_emps_lbx;
    logs{thread_index}.ubx = nlid_emps_ubx;
    logs{thread_index}.x0_est_v = x0_est_v;
    [~,~,logs{thread_index}.rel_errors_training] = nlid_emps_sim_plant(x0_est_v, u_data_v{initial_params_index}, param_est_comb, y_data_v{initial_params_index}); %against the data we used to train
    diary nlid_emps_sweep_logs.txt; disp(['hybrid_group_size = ' num2str(hybrid_group_size) ...
        ' initial_params_index = ' num2str(initial_params_index) ... 
        ' time_toc = ' num2str(logs{thread_index}.time_toc) ... 
        ' rmserror = ' num2str(rmserror) ...
        ' param_est_comb = ' num2str(param_est_comb') ...
        ' x0_est_v = ' num2str(x0_est_v.') ...
        ]); 
        diary off
end
save('nlid_emps_sweep_logs.mat','logs','u_data_v', 'y_data_v','x0_v','X_sim_H','param_guess_v')

%% Plot: number of cases with rel. error lower than 0.1%
to_plot=[]
for log_item = logs
    to_plot(log2(log_item{1}.hybrid_group_size)+1, log_item{1}.initial_params_index) = log_item{1}.rel_errors_training.plant_qm;
end
colororder(repmat(linspace(0,0.8,7),3,1)')
clf
bar(sum(to_plot.'<0.1))
set(gca,'LineWidth',2)
ylabel("Number of successful cases")
xticks(1:18)
xticklabels({"MS | D=1","D=2","D=4","D=8","D=16","D=32","D=64","D=128","D=256","D=512","D=1024","D=2048","D=4096","D=8192","SS | D=10000","SS | nlgreyest","PEM+SS","PEM"})
xtickangle(90)
set(gca,'FontSize',24)
exportgraphics19 paperplot_num_success

%% Plot: running time for successful cases
to_plot=[];
for log_item = logs
    to_plot(log2(log_item{1}.hybrid_group_size)+1, log_item{1}.initial_params_index) = iif(log_item{1}.rel_errors_training.plant_qm<0.1, log_item{1}.time_toc, nan);
end
colororder(repmat(linspace(0,0.8,7),3,1)')
clf
hlc=boxplot(to_plot.');
xticks(1:18)
xticklabels({"MS | D=1","D=2","D=4","D=8","D=16","D=32","D=64","D=128","D=256","D=512","D=1024","D=2048","D=4096","D=8192","SS | D=10000","SS | nlgreyest","PEM+SS","PEM"})
xtickangle(90)
ylabel('time (seconds)')
set(gca,'FontSize',24)
set(gca,'LineWidth',2)
for ih=1:7, set(hlc(ih,:),'LineWidth',2);, end
exportgraphics19 paperplot_time_toc

%% Plot: relative error for successful cases
to_plot = []
for log_item = logs
    to_plot(log2(log_item{1}.hybrid_group_size)+1, log_item{1}.initial_params_index) = iif(log_item{1}.rel_errors_training.plant_qm<0.1, log_item{1}.rel_errors_training.plant_qm, nan);
end
colororder(repmat(linspace(0,0.8,7),3,1)')
clf
hlc=boxplot(to_plot.');
xticklabels({"MS | D=1","D=2","D=4","D=8","D=16","D=32","D=64","D=128","D=256","D=512","D=1024","D=2048","D=4096","D=8192","SS | D=10000","SS | nlgreyest","PEM+SS","PEM"})
xtickangle(90)
set(gca,'FontSize',24)
set(gca,'LineWidth',2)
ylabel('relative error (%)')
for ih=1:7, set(hlc(ih,:),'LineWidth',2); end
exportgraphics19 paperplot_rel_error
