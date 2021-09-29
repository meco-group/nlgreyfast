%% Identify a mass-spring-damper
% In this example, we will simulate the output of a mass-spring-damper system with known parameters (m, b and k),
% and will use different formulations of the optimal control problem (single shooting, multiple shooting, hybrid) 
% to identify those parameters, fitting its output to the previously simulated data.
% Variable name postfixes _v and _h correspond to vertical and horizontal vectors, 
% _V corresponds to a "tall" matrix (the data records are in rows), _H corresponds to a "long" matrix (the data records are in columns).
clear all
close all

%% Settings
N = 10000;  % Number of samples that we will generate in simulation
fs = 1000; % Sampling frequency [Hz]
% To make the RK4 integration accurate, we will take intermediate steps in between the samples while calculating 
% the output of the system for a given input:
N_steps_per_sample = 10; 
% These are the known parameters [m;b;k] of the mass-spring-damper:
param_truth_v = [0.1;0.2;1];
% The solver will start from this guess of the parameters:
param_guess_v = [0.08;0.08;0.5];
% We can rescale the parameters to the same range (if needed) so that the solver will not fail due to numerical 
% problems. By default we apply no scaling: [1;1;1]
scale_v = [1;1;1];

%% Modeling
% At this point we will input the ODE for the mass-spring-damper:

% States:
y  = casadi.MX.sym('y'); 
dy = casadi.MX.sym('dy');
states_v = [y;dy];

% Inputs to the system:
u  = casadi.MX.sym('u'); %u is the input
controls_v = u;

% The parameters to estimate:
m = casadi.MX.sym('m');
b = casadi.MX.sym('b');
k = casadi.MX.sym('k');
params_v = [m;b;k];

% The ODEs are the following: 
% \dot y = dy
% \dot dy = (-k*dy-b*y+u)/m
% The functions in this toolbox expect the ODE to be specified as a CasADi function that returns the right hand side
% of the ODE, depending on the `states`, `controls` as the inputs, and `params` as the parameters to estimate.
rhs = [dy; (-k*dy-b*y+u)/m];
ode = casadi.Function('ode',{states_v,controls_v,params_v},{rhs}); % Form an ODE function

%% Simulate the system
% In this step we simulate the mass-spring-damper, given a random signal as input/excitation:
u_data_v = 0.1*rand(N,1); 
% ...starting from the following states:
x0_v = [0;0];
% Now we formulate a CasADi function that can perform the simulation (knowing the ones above):
[plant_next,all_samples] = nlsim(ode, fs, N-1, N_steps_per_sample); 
% Actually compute the simulation for x0 and u_data, for the known real parameters:
X_measured_H = full([x0_v all_samples(x0_v.', u_data_v(1:end-1).', repmat(param_truth_v,1,N-1))]); 
% Only some of the states will be measured / included in the output, their indexes are listed in the following array:
y_states_v = [1];
y_data_v = X_measured_H(y_states_v,:).'; %now actually select those states

% We may add some noise to this simulated output if we uncomment the line below. When noise is absent, the fit will be perfect.
%y_data_v = y_data_v + 0.001 * randn(size(y_data_v))

%% Single shooting / multiple shooting / hybrid shooting
% N_group_size = N; corresponds to single shooting formulation
% N_group_size = 1; corresponds to multiple shooting formulation
% N_group_size = 4; and any other integer value between 1 and N correspond to the hybrid shooting / partially 
%   constrained multiple shooting formulation
N_group_size = 4; 

[param_est_comb_v, x0_est_comb_v, diag_data_comb] = nlidcomb( ...
    plant_next, ...     % The system to be identified
    fs, ...             % The sampling rate (important for 'guess_mechatronic' only)
    u_data_v, ...       % The measured input
    y_data_v, ...       % The measured output
    y_states_v, ...     % Which state corresponds to the output (only applicable if y_mapping is empty), thus y=x(y_states) is assumed for y_mapping
    [], ...             % y_mapping to map from states to output (the y=G(x,u) part of the NLSS model) 
    [0;0;0], ...        % lower bounds for params to estimate 
    param_guess_v, ...  % initial guess for params to estimate
    [inf;inf;inf], ...  % upper bounds for params to estimate
    [], ...             % lower bounds for the elements to be estimated in x0
    x0_v, ...           % initial guess for x0
    [], ...             % upper bounds for the elements to be estimated in x0
    [false;false], ...  % Which of the initial states x0 are to be included in the estimation. Note that it can make sense to estimate them because even if we know them from measurement, they can be noisy.)
    scale_v, ...        % Explained above.
    N_group_size, ...   % >>> Probably the most important parameter, controls the recursion depth / group size. See explanation above.
    [], ...             % lower bounds for how the states evolve (if X_guess_H is 'double' then it should be the same size as X_guess_H or [] to be automatically initialized) 
    'guess_mechatronic', ... % X_guess_H: initial guess for how the states evolve, see notes on 'guess_mechatronic' later
    [], ...             % upper bounds for how the states evolve (if X_guess_H is 'double' then it should be the same size as X_guess_H or [] to be automatically initialized) 
    [], ...             % additional constraints lower bound
    [], ...             % expression for additional constraints
    [], ...             % additional constraints upper bound
    'all', ...          % Which of the gap constraints should be taken into account for MS/hybrid method. 
    '');                % If we provide '-march=native -Ofast' here, then JIT will be switched on: CasADi will attempt
                        % to precompile the expression with GCC + ccache to make it faster to evaluate.
                        % However, CasADi will not be able to do thread-based parallelisation of multiple shooting.

% ## Notes on 'guess_mechatronic'
% In the MS and hybrid method, we have to give an initial guess on how the states evolve with time (the solver will later refine it). 
%     In a mechatronic system, where the ODE is like
%         \dot q = dq
%         \dot dq = something
%     if we measure all the q (as in the case of MSD), we can just use finite differences as a guess.
%     This specific case is denoted with 'guess_mechatronic'.
%     There are also other possibilities to pass as this parameter, e.g. you can provide a matrix with the values for the guesses,
%     or provide 'guess_exact' for a system which is like \dot dq = something.

disp(['||estimated parameters - known parameters||₂ = ' num2str(norm(param_est_comb_v-param_truth_v,'inf'))])
