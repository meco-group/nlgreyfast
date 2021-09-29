function estmodel = nlgreyfast(data, model, fopts)
    % Tool for nonlinear grey-box system identification based on CasADi
    if strcmp(class(data),'iddata')==1
        u_data_V = data.InputData;
        y_data_V = data.OutputData;
        assert(data.Ts~=0,'only discrete time systems are supported')
        fs = 1/data.Ts;
    else
        error('argument ''data'' should be iddata')
    end
    if strcmp(class(model),'casadi.Function')
        plant_next = model;
    else
        error('argument ''model'' should be casadi.Function')
    end

    % in fopts, at minimum you should define: param_guess_v, x0_v, N_group_size, X_guess_H 
    if ~isfield(fopts, 'param_guess_v'), error('fopts.param_guess_v not defined'), end
    if ~isfield(fopts, 'x0_v'), error('fopts.x0_v not defined'), end
    if ~isfield(fopts, 'param_lb_v'), fopts.param_lb_v = -inf(size(fopts.param_guess_v)); end
    if ~isfield(fopts, 'param_ub_v'), fopts.param_ub_v =  inf(size(fopts.param_guess_v)); end
    if ~isfield(fopts, 'y_mapping'), fopts.y_mapping = []; end
    if ~isfield(fopts, 'lbx_x0_free_v'), fopts.lbx_x0_free_v = []; end
    if ~isfield(fopts, 'ubx_x0_free_v'), fopts.ubx_x0_free_v = []; end
    if ~isfield(fopts, 'x0_free_v'), fopts.x0_free_v = true(size(fopts.x0_v)); end
    if ~isfield(fopts, 'scale_v'), fopts.scale_v = ones(size(fopts.param_guess_v)); end
    if ~isfield(fopts, 'N_group_size'), error('fopts.N_group_size not defined'), end
    if ~isfield(fopts, 'lbx_X_guess_H'), fopts.lbx_X_guess_H = []; end
    if ~isfield(fopts, 'X_guess_H'), error('fopts.X_guess_H not defined'), end
    if ~isfield(fopts, 'ubx_X_guess_H'), fopts.ubx_X_guess_H = []; end
    if ~isfield(fopts, 'additional_constraints_lbg_v'), fopts.additional_constraints_lbg_v = []; end
    if ~isfield(fopts, 'additional_constraints_g'), fopts.additional_constraints_g = []; end
    if ~isfield(fopts, 'additional_constraints_ubg_v'), fopts.additional_constraints_ubg_v = []; end
    if ~isfield(fopts, 'enable_gaps_v'), fopts.enable_gaps_v = 'all'; end

    estmodel = idnlgreyfast;
    [estmodel.param_est_v, estmodel.x0_est_v, estmodel.diag_data ] = nlidcomb(plant_next, fs, u_data_V, y_data_V, ...
        fopts.y_states_v, fopts.y_mapping, fopts.param_lb_v, fopts.param_guess_v, fopts.param_ub_v, ...
        fopts.lbx_x0_free_v, fopts.x0_v, fopts.ubx_x0_free_v, fopts.x0_free_v, fopts.scale_v, fopts.N_group_size, ...
        fopts.lbx_X_guess_H, fopts.X_guess_H, fopts.ubx_X_guess_H, fopts.additional_constraints_lbg_v, fopts.additional_constraints_g, ...
        fopts.additional_constraints_ubg_v, fopts.enable_gaps_v, fopts);

end