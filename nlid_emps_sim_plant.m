function [qm_sim_v X_sim_H rel_errors] = nlid_emps_sim_plant(x0_v, vir_v, param_sim_v, qm_measured_v)
    % simulate EMPS for given input
    % test this function: 
    %   [qm_sim_v, X_sim_V] = nlid_emps_sim_plant(x0_plant, vir,  param_truth);
    %   dlot(X_sim_V);

    assert_v(x0_v)
    assert_v(param_sim_v)
    assert_v(vir_v)
    [plant_next, plant_ddq_out, system_next, system_vir_out, system_ddq_out, controller_next, controller_vir_out, ode, plant_states, plant_controls, system_states, system_controls, params_sym, system_pulse_next, system_pulse_states, system_pulse_controls, system_pulse_vir_out, system_pulse_ddq_out, tec, fs, system_noise_next] = nlid_emps_models;
    N=size(vir_v,1);
    plant_next_all_samples = plant_next.mapaccum('plant_next_all_samples', N-1);
    X_sim_H = full([x0_v plant_next_all_samples(x0_v.',vir_v(1:end-1).',repmat(param_sim_v,1,N-1))]);
    X_sim_V = X_sim_H.';
    qm_sim_v = X_sim_V(:,1);
    rel_errors = struct;
    if nargin >= 4
        assert_v(qm_measured_v)
        rel_errors.plant_qm = 100*norm(qm_measured_v - qm_sim_v)/norm(qm_measured_v);
    end