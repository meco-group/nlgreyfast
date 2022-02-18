function estmodel = nlgreyfast(data, model, fopts)
    % Tool for nonlinear grey-box system identification based on CasADi.
    % nlgreyfast is a wrapper around the function nlidcomb.
    % Output:
    % - *estmodel*: an idnlgreyfast object (contains estimated parameters, initial states, and debug data; currently not compatible with the idnlgrey object of MATLAB, i.e. you cannot use compare() and other functions on it from MathWorks System Identification Toolbox).
    % Input parameters:
    % - *data*: an *iddata* object from MathWorks System Identification Toolbox. (Only its InputData, OutputData and Ts are parsed.)
    % - *model*: a CasADi function that tells the next state based on the current state
    % - *fopts*: a struct with various options, at minimum you should define: param_guess_v, x0_v, N_group_size, X_guess_H 
    %   - *param_guess_v* (required): the initial parameter values for all free parameters
    %   - *x0_v* (required): the initial state values for all initial states, including free and fixed ones
    %   - *y_states_v* (required): integer numbers to specify which of the states can be measured directly (i.e. these will be the outputs). If you want to specify a state output equation, like y=g(x,u,theta) in the NL state-space model, then you need to use *y_mapping*. *nlidcomb* can use that some of the states are measurable directly (e.g. for generating *X_guess*, or for PEM-like formulations, see  *fix_states_to_X_guess_H_v*), that's also a benefit of *y_states_v*. See also *X_guess_H* for how this is used.
    %   - *param_lb_v*: the lower bounds for all free parameters.
    %   - *param_ub_v*: the upper bounds for all free parameters.
    %   - *y_mapping*: _mapping to map from states to output (the y=G(x,u) part of the NLSS model).
    %   - *lbx_x0_free_v*: the lower bounds for all free initial states.
    %   - *ubx_x0_free_v*: the upper bounds for all free initial states.
    %   - *x0_free_v*: true or false values to define which of the elements in x0 are free (true) or fixed (false).
    %   - *scale_v*: allows to rescale the parameters to the same range (if needed) so that the solver will not fail due to numerical problems. If not defined, we apply no scaling: [1;1; ... ;1]
    %   - *N_group_size* (required): the subrecord length / recursion depth of the (partially constrained) multiple shooting. Value 1 means multiple shooting, if value equals number of input samples, then we have the single shooting formulation, in between we have the partially constrained multiple shooting formulations.
    %   - *lbx_X_guess_H*: lower bounds for how the states evolve (if X_guess_H is 'double' then it should be the same size as X_guess_H or [] to be automatically initialized).
    %   - *X_guess_H* (required): initial guess for how the states evolve. Possible values:
    %      - *guess_mechatronic*: in the state vector, first we have some states that we can directly measure, like positions. These measurements are available in *y_data_V* and will be used for the states marked true in *y_states_v*. The states afterwards are initialized using finite differences (diff()) from these states, so that they are the velocities corresponding to the positions, using the same example. If there are any remaining states, they are initialized as 0.
    %      - *guess_exact*: similar to *guess_mechatronic*, but no states are initialized with finite differences, and no remaining states are assumed.
    %      - a matrix of floating point numbers that should have the same length (dim 2) as the input and output signals specified, and the same height (dim 1) as the number of states: we can directly specify the guess for how the states evolve.
    %   - *ubx_X_guess_H*: upper bounds for how the states evolve (if X_guess_H is 'double' then it should be the same size as X_guess_H or [] to be automatically initialized) 
    %   - *additional_constraints_lbg_v*: additional constraints lower bound
    %   - *additional_constraints_g*: additional constraints as CasADi function. It will be treated as: additional_constraints_lbg_v <= additional_constraints_g(X_sim_H,u_data_H,parameters) <= additional_constraints_ubg_v 
    %   - *additional_constraints_ubg_v*: additional constraints upper bound
    %   - *enable_gaps_v*: allows you to enable or disable gap constraints (default: 'all', or takes a vector of true or false values; [] will disable all gap constraints).
    %   - *fit_against_X_guess_H_v*: integer numbers to specify which of the states shall be fixed to *X_guess_H* at the beginning of subrecords. These will not have a corresponding decision variable in any formulation. It is useful for specifying PEM-like methods.
    %   - *fix_states_to_X_guess_H_v*: integer numbers to specify which of the states shall be fit against *X_guess_H*. These are additional terms in the objective. It is useful for specifying PEM-like methods.


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

    % in fopts, at minimum you should define: param_guess_v, x0_v, N_group_size, X_guess_H, y_states_v 
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
