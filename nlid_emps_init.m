% Initialize workspace with EMPS related variables

clear all
close all
clc
delete(gcp('nocreate'))

%%
load('nlid_emps_sim_data_training.mat')
M1  = 95.1089;
Fv1 = 203.5034;
Fc1 = 20.3935;
OF1 = -3.1648;

%% Compute relative errors on original data
rel_err_q   = 100*norm(q_f(n_border:N_pts) - q_s(n_border:N_pts))/norm(q_f(n_border:N_pts));
rel_err_dq  = 100*norm(dq_f(n_border:N_pts) - dq_s(n_border:N_pts))/norm(dq_f(n_border:N_pts));
rel_err_ddq = 100*norm(ddq_f(n_border:N_pts) - ddq_s(n_border:N_pts))/norm(ddq_f(n_border:N_pts));
rel_err_u   = 100*norm(Force1(n_border:N_pts) - Force1_s(n_border:N_pts))/norm(Force1(n_border:N_pts));
rel_err_qm   = 100*norm(qm - q_s)/norm(qm);

%% Build all the models with nlid_emps_models function
param_truth = [M1; Fv1; Fc1; OF1];
[plant_next, plant_ddq_out, system_next, system_vir_out, system_ddq_out, controller_next, controller_vir_out, ode, plant_states, plant_controls, system_states, system_controls, params_sym, system_pulse_next, system_pulse_states, system_pulse_controls, system_pulse_vir_out, system_pulse_ddq_out, tec, fs, system_noise_next, system_newcontroller_next, system_newcontroller_states, system_newcontroller_controls, system_newcontroller_vir_out, system_noise_vir_out, plant_substm_next] = nlid_emps_models;
x0_system = [q_f(1);dq_f(1);0;0;0];
x0_plant=[qm(1);dq_f(1)];
x0_controller = [0;0;0];

%% Take the true parameters to simulate 
param_sim = param_truth;

%% Simulate everything using CasADi
%nlid_emps_sim(param_truth)

%% Simulate and plot everything using CasADi
%close all, nlid_emps_sim(param_truth,'plot',true)

%% Initialization for identification -- wild bounds 
ubx = [inf(3,1);0];
lbx = [0;0;0;-inf];
param_guess=1e2*[1;1;1;-1];
y_states = [1];
y_mapping = casadi.Function('y_mapping',{system_states,system_controls,params_sym},{veccat(system_states{1}, system_vir_out(system_states,system_controls,params_sym))});
scale = ones(size(param_truth));
additional_lbg = repmat(-10,length(qm),1);
additional_ubg = repmat(10, length(qm),1);
additional_g = system_vir_out.map(length(qm));
%X_guess = full(system_next_all_samples(x0_system,qg,repmat(param_guess,1,N)));