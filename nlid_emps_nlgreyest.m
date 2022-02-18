function [param_est_nlgreyest, x0_est, nlgr] = nlid_emps_nlgreyest(EvalFun, param_guess_v, x0_v, fopts)
    % Estimate the parameters and initial state of the EMPS using nlgreyest (from the MathWorks System Identification Toolbox).
    % Inputs: 
    % - *EvalFun*: the name of a MATLAB function defined in an M-file, that contains the state and output equations of the model (see *nlgreyest* documentation).
    % - *param_guess_v*: initial parameters for the estimation
    % - *x0_v*: initial states for the simulation of the states (will not be decision variables)
    % - *fopts*: struct with additional options:
    %    - *param_lb_v* and *param_ub_v*: lower and upper bounds for the parameters to be estimated
    %    - *x0_free_v*: vector of true and false values: which of the initial states should be estimated (true) or not (false)
    %    - *u_data_v* and *y_data_v* allows to set custom input and output data for parameter estimation, otherwise they'll be taken from *nlid_emps_sim_data_training.mat*. 
    % 
    % Command to test it: 
    %   nlid_emps_init
    %   fopts = struct; fopts.x0_free_v = [false;true];
    %   [param_est_nlgreyest, x0_est_nlgreyest, nlgr_nlgreyest] = nlid_emps_nlgreyest('emps_plant_atan_m_mex', param_guess, x0_plant, fopts)
    
    if ~isfield(fopts, 'param_lb_v'), fopts.param_lb_v = -inf(4,1); end
    if ~isfield(fopts, 'param_ub_v'), fopts.param_ub_v = inf(4,1); end
    if ~isfield(fopts, 'x0_free_v'), fopts.x0_free_v = false(2,1); end

    load('nlid_emps_sim_data_training');
    if ~isfield(fopts, 'u_data_v'), fopts.u_data_v = vir; end
    if ~isfield(fopts, 'y_data_v'), fopts.y_data_v = qm; end
    
    pas = 9.999839611217580e-04;

    iddata_emps_plant = iddata(fopts.y_data_v, fopts.u_data_v, pas, 'Name', 'EMPS');
    iddata_emps_plant.InputName = 'vir'; %if we have multiple, then we can use {}
    iddata_emps_plant.InputUnit = 'V';
    iddata_emps_plant.OutputName = 'qm';
    iddata_emps_plant.Tstart = 0;
    iddata_emps_plant.TimeUnit = 's';

    Parameters = [];
    Parameters(1).Value = param_guess_v(1);
    Parameters(1).Name = 'M1'; 
    Parameters(1).Minimum = fopts.param_lb_v(1);
    Parameters(1).Maximum = fopts.param_ub_v(1);
    Parameters(1).Unit = ''; 
    Parameters(1).Fixed = false;
    Parameters(2).Value = param_guess_v(2);
    Parameters(2).Name = 'Fv1'; 
    Parameters(2).Minimum = fopts.param_lb_v(2); 
    Parameters(2).Maximum = fopts.param_ub_v(2); 
    Parameters(2).Unit = ''; 
    Parameters(2).Fixed = false; 
    Parameters(3).Value = param_guess_v(3);
    Parameters(3).Name = 'Fc1'; 
    Parameters(3).Minimum = fopts.param_lb_v(3); 
    Parameters(3).Maximum = fopts.param_ub_v(3); 
    Parameters(3).Unit = ''; 
    Parameters(3).Fixed = false; 
    Parameters(4).Value = param_guess_v(4);
    Parameters(4).Name = 'OF1'; 
    Parameters(4).Minimum = fopts.param_lb_v(4); 
    Parameters(4).Maximum = fopts.param_ub_v(4); 
    Parameters(4).Unit = ''; 
    Parameters(4).Fixed = false;
    Ts = 0; %Continuous-time system (EvalFun interpreted as ODE: dx=f(x))
    EvalArgument = {};

    Order = [1 1 2]; % Model orders [ny nu nx].
    InitialStates = x0_v; % Initial initial states.
    nlgr = idnlgrey(EvalFun, Order, Parameters, InitialStates, Ts, 'Name', 'EMPS','FileArgument', EvalArgument);
    set(nlgr, 'InputName', 'vir', 'InputUnit', 'V', 'OutputName', 'qm', 'OutputUnit', 'm', 'TimeUnit', 's');
    %nlgr.SimulationOptions.AbsTol = 1e-6;
    %nlgr.SimulationOptions.RelTol = 1e-5;
    %compare(iddata_emps_plant, nlgr);

    opt = nlgreyestOptions('Display', 'full','SearchMethod','lsqnonlin');
    nlgr = setinit(nlgr, 'Fixed', {~fopts.x0_free_v(1) ~fopts.x0_free_v(2)}); %here the meaning is fixed, not free, hence the negation
    nlgr = nlgreyest(iddata_emps_plant, nlgr, opt);
    nlgr.Report
    fprintf('\n\nThe search termination condition:\n')
    nlgr.Report.Termination

    param_est_nlgreyest = [];
    param_est_nlgreyest(1,1) = nlgr.Parameters(1).Value;
    param_est_nlgreyest(2,1) = nlgr.Parameters(2).Value;
    param_est_nlgreyest(3,1) = nlgr.Parameters(3).Value;
    param_est_nlgreyest(4,1) = nlgr.Parameters(4).Value;
    
    x0_est = [nlgr.InitialStates(1).Value; nlgr.InitialStates(2).Value];
end
