% Example for identification of the EMPS system.
% The EMPS consists of a DC motor connected to a ballscrew drive with a load, the position of 
% which is measured by an encoder. The model is a SISO, nonlinear state-space model with two states.
% - Input: control voltage of the motor.
% - Output: position of the load.
% - States: position and velocity of the load.
% The four parameters to estimate: load inertia, viscous and Coluomb friction coefficients, measurement offset.
% See our paper for a more detailed description of the system and a figure.
% The _v, _h, _V, _H postfixes contain information about the size of the matrix, see `nlidcomb.h`.

clear all
close all
clc

% These were the parameters acquired with the IDIM-LS method with the scripts of the original paper on the EMPS:
param_idimls_v = [95.1089; 203.5034; 20.3935; -3.1648]; 

% Load data
load('nlid_emps_sim_data_training.mat')
model = nlid_emps_models;

%% Identify EMPS using PCMS.
% The subrecord length should be defined here:
D = 8;             % multiple shooting
%D = size(qm,1);   % single shooting
%partially constrained multiple shooting if 1 < D < number of samples

fopts = struct;
fopts.scale_v = ones(4,1); %We expect all parameters to be in the same range.
fopts.param_guess_v = [100;100;100;-100]; %This is the initial guess for the parameters to estimate.
fopts.x0_v = [qm(1);dq_f(1)]; %This is the initial state of the estimation.
fopts.x0_free_v = [false;false]; %Shall we include the initial states into the set of decision variables? If so, which ones?
fopts.y_states_v = [1]; %The first state is an output.
fopts.param_lb_v = [0;0;0;-inf]; %Lower bound for the parameters to estimate.
fopts.param_ub_v = [inf;inf;inf;0]; %Upper bound for the parameters to estimate.
fopts.X_guess_H = 'guess_mechatronic'; %Initializing the decision variables corresponding to how the states evolve.
    % In this case we have the measurements of the first state (as denoted in `y_states_v`).
    % The second state will be inferred from the first one through finite differences.
    % The measured / inferred values will be used to initialize the corresponding decision variables.
fopts.N_group_size = D; %The subrecord length.

Ts = 9.999839611217580e-04; % The sample time.
data = iddata(qm, vir, Ts); % qm: output, vir: input

estmodel_pcms = nlgreyfast(data, model, fopts); % Do the estimation

disp('PCMS -- estimated parameters:')
disp(estmodel_pcms.param_est_v)
disp('PCMS -- estimated initial states (if fixed, then the same as the initial):')
disp(estmodel_pcms.x0_est_v)


%% Identify EMPS using PEM+SS
fopts.fit_against_X_guess_H_v = [1; 2];
fopts.fix_states_to_X_guess_H_v = [1; 2];
fopts.N_group_size = 1;
estmodel_pem = nlgreyfast(data, model, fopts); 
fopts.x0_v = estmodel_pem.x0_est_v;
fopts.param_guess_v = estmodel_pem.param_est_v;
fopts.fit_against_X_guess_H_v = [];
fopts.fix_states_to_X_guess_H_v = [];
fopts.N_group_size = size(qm,1);
estmodel_pemss = nlgreyfast(data, model, fopts); 

disp('PEM+SS -- estimated parameters:')
disp(estmodel_pemss.param_est_v)
disp('PEM+SS -- estimated initial states (if fixed, then the same as the initial):')
disp(estmodel_pemss.x0_est_v)

%% Plots of measured vs. simulated system output
[qm_sim_idimls_v, ~, rel_errors_idimls] = nlid_emps_sim_plant(fopts.x0_v, vir, param_idimls_v, qm);
[qm_sim_pcms_v, ~, rel_errors_pcms] = nlid_emps_sim_plant(estmodel_pcms.x0_est_v, vir, estmodel_pcms.param_est_v, qm);
[qm_sim_pemss_v, ~, rel_errors_pemss] = nlid_emps_sim_plant(estmodel_pemss.x0_est_v, vir, estmodel_pemss.param_est_v, qm);
clf
hold on
plot(qm)
plot(qm_sim_idimls_v)
plot(qm_sim_pcms_v)
plot(qm_sim_pemss_v)
title('measured vs. simulated system output')
xlabel('#N samples')
ylabel('position')
legend('measured output', ['IDIM-LS|error: ' num2str(rel_errors_idimls.plant_qm) '%'], ['PCMS|error: ' num2str(rel_errors_pcms.plant_qm) '%'], ['PEM+SS|error: ' num2str(rel_errors_pemss.plant_qm) '%'])
hold off