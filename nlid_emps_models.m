function [plant_next, plant_ddq_out, system_next, system_vir_out, system_ddq_out, controller_next, controller_vir_out, ode, plant_states, plant_controls, system_states, system_controls, params_sym, system_pulse_next, system_pulse_states, system_pulse_controls, system_pulse_vir_out, system_pulse_ddq_out, tec, fs, system_noise_next, system_newcontroller_next, system_newcontroller_states, system_newcontroller_controls, system_newcontroller_vir_out, system_noise_vir_out, plant_substm_next] = nlid_emps_models
    %% Build all the ODEs for the EMPS in CasADi. There are several versions of the ODEs corresponding to different experiments.
    import casadi.*

    %internal parameters of this function:
    N_steps_per_sample = 10; %interpolation steps

    % data we use from the EMPS model files:
    gtau = 35.150651882485470;
    filterq1 = [0.5000 0.5000];
    tec = 9.999839611217580e-04;
    %tec=1e-3;
    fs = 1/tec;
    kv = 2.434500000000000e+02;
    kp = 1.601800000000000e+02;

    %% [PLANT]
    %states:
    plant_q  = MX.sym('plant_q',1);
    plant_dq  = MX.sym('plant_dq',1);
    plant_states = veccat(plant_q,plant_dq);

    %controls:
    plant_u_vir  = MX.sym('plant_u_vir',1);
    plant_controls = plant_u_vir;

    %parameters:
    M1_sym = MX.sym('M1_sym',1);
    Fv1_sym = MX.sym('Fv1_sym',1);
    Fc1_sym = MX.sym('Fc1_sym',1);
    OF1_sym = MX.sym('OF1_sym',1);
    params_sym = veccat(M1_sym, Fv1_sym, Fc1_sym, OF1_sym);

    %friction formula used to build plant_ddq:
    fric_simulink_formula_sym = @(x)sign(x)*(Fv1_sym*abs(x)+Fc1_sym); %this one is not differentiable due to the sign()
    fric_continuous_formula_sym = @(x)Fc1_sym*tanh(10000/Fc1_sym*x)+Fv1_sym*x;
    friction = Function('friction',{plant_dq,Fc1_sym,Fv1_sym},{fric_continuous_formula_sym(plant_dq)}); %can choose between two formulas above

    %plant ODE
    plant_ddq = (1/M1_sym)*(gtau*plant_controls - OF1_sym - friction(plant_dq,Fc1_sym,Fv1_sym));
    rhs = veccat(plant_dq, plant_ddq);
    ode = Function('ode',{plant_states,plant_controls,params_sym},{rhs});
    plant_next = nlsim(ode, fs, [], N_steps_per_sample); %next state of plant
    plant_ddq_out = Function('plant_ddq_out',{plant_states,plant_controls,params_sym},{plant_ddq}); %output of plant (based on current sample)

%% [PLANT+plant_substm_next] 
    %plant ODE
    M_subst = M1_sym-(1/(M1_sym-100)^200);
    plant_ddq_substm = (1/M_subst)*(gtau*plant_controls - OF1_sym - friction(plant_dq,Fc1_sym,Fv1_sym));
    rhs_substm = veccat(plant_dq, plant_ddq_substm);
    ode_substm = Function('ode',{plant_states,plant_controls,params_sym},{rhs_substm});
    plant_substm_next = nlsim(ode_substm, fs, [], N_steps_per_sample); %next state of plant


    %% [SYSTEM + PULSE]
    input_qg = MX.sym('input_qg');
    pulse = MX.sym('pulse');

    %plant input
    unit_delay_state = MX.sym('unit_delay_state');
    system_pulse_controls = veccat(input_qg,pulse);
    clamping_value=10;
    clamping_formula = @(x)min(max(-clamping_value,x),clamping_value);
    sat_vir_out = clamping_formula(unit_delay_state)+pulse; %we have clamping on vir for more accurate simulation... this formula will not be used in optimization though
    %sat_vir_out = unit_delay_state; %disabling clamping on vir

    %plant output
    q1_vect = plant_next(plant_states,sat_vir_out,params_sym);

    %filter block
    filter_state = MX.sym('filter_state');
    filter_out = filterq1(1)*q1_vect(1)+filterq1(2)*filter_state;
    filter_state_next = q1_vect(1);

    %derivative block
    derivative_state = MX.sym('derivative_state');
    derivative_out = (filter_out-derivative_state)/tec;
    derivative_state_next = filter_out;

    %unit delay
    unit_delay_state_next = kv * (kp*(input_qg-q1_vect(1))-derivative_out);

    %system formulas
    system_pulse_states = veccat(plant_states, unit_delay_state, derivative_state, filter_state);
    system_pulse_next = Function('system_pulse_next',{system_pulse_states,system_pulse_controls,params_sym}, {veccat( ...
        q1_vect, unit_delay_state_next, derivative_state_next, filter_state_next ...
        )},{'x','u','m'},{'x_next'}); %next state of system
    system_pulse_vir_out = Function('system_pulse_vir_out', {system_pulse_states,system_pulse_controls,params_sym}, {sat_vir_out}); %system output: control voltage of plant
    system_pulse_ddq_out = Function('system_pulse_ddq_out',{system_pulse_states,system_pulse_controls,params_sym},{plant_ddq_out(system_pulse_states(1:2),system_pulse_vir_out(system_pulse_states,system_pulse_controls,params_sym),params_sym)}); %system output: acceleration output of the plant
    

    %% [SYSTEM]
    input_qg = MX.sym('input_qg');

    %plant input
    unit_delay_state = MX.sym('unit_delay_state');
    system_controls = input_qg;
    %sat_vir_out = clamping_formula(unit_delay_state);
    sat_vir_out = unit_delay_state; %disabling clamping on vir

    %plant output
    q1_vect = plant_next(plant_states,sat_vir_out,params_sym);

    %filter block
    filter_state = MX.sym('filter_state');
    filter_out = filterq1(1)*q1_vect(1)+filterq1(2)*filter_state;
    filter_state_next = q1_vect(1);

    %derivative block
    derivative_state = MX.sym('derivative_state');
    derivative_out = (filter_out-derivative_state)/tec;
    derivative_state_next = filter_out;

    %unit delay
    unit_delay_state_next = kv * (kp*(input_qg-q1_vect(1))-derivative_out);

    %system formulas
    system_states = veccat(plant_states, unit_delay_state, derivative_state, filter_state);
    system_next = Function('system_next',{system_states,system_controls,params_sym}, {veccat( ...
        q1_vect, unit_delay_state_next, derivative_state_next, filter_state_next ...
        )},{'x','u','m'},{'x_next'}); %next state of system
    system_vir_out = Function('system_vir_out', {system_states,system_controls,params_sym}, {sat_vir_out}); %system output: control voltage of plant
    system_ddq_out = Function('system_ddq_out',{system_states,system_controls,params_sym},{plant_ddq_out(system_states(1:2),system_vir_out(system_states,system_controls,params_sym),params_sym)}); %system output: acceleration output of the plant

    
    %% [CONTROLLER] (for simulating the controller while substituting the plant)
    input_q = MX.sym('input_q');
    input_qg = MX.sym('input_qg');
    controller_controls = veccat(input_qg, input_q);

    % we alter the formulas for filter, derivative and unit delay blocks so that it's without [PLANT]
    filter_out = filterq1(1)*input_q+filterq1(2)*filter_state;
    filter_state_next = input_q;
    derivative_out = (filter_out-derivative_state)/tec;
    derivative_state_next = filter_out;
    unit_delay_state_next = kv * (kp*(input_qg-input_q)-derivative_out);

    %controller formulas
    controller_states = veccat(unit_delay_state, derivative_state, filter_state);
    controller_next = Function('controller_next',{controller_states,controller_controls}, {veccat( ...
        unit_delay_state_next, derivative_state_next, filter_state_next ...
        )},{'x','u'},{'x_next'});
    controller_vir_out = Function('controller_vir_out', {controller_states,controller_controls}, {sat_vir_out}); %simulate the control voltage (plant input)


    %% [SYSTEM+noise]
    input_noise = MX.sym('input_noise');
    input_qg = MX.sym('input_qg');

    %plant input    
    unit_delay_state = MX.sym('unit_delay_state');
    system_noise_controls = veccat(input_qg, input_noise);
    %sat_vir_out = clamping_formula(unit_delay_state);
    sat_vir_out = unit_delay_state; %disabling clamping on vir

    %plant output
    q1_vect = plant_next(plant_states,sat_vir_out,params_sym)+input_noise;

    %filter block
    filter_state = MX.sym('filter_state');
    filter_out = filterq1(1)*q1_vect(1)+filterq1(2)*filter_state;
    filter_state_next = q1_vect(1);

    %derivative block
    derivative_state = MX.sym('derivative_state');
    derivative_out = (filter_out-derivative_state)/tec;
    derivative_state_next = filter_out;

    %unit delay
    unit_delay_state_next = kv * (kp*(input_qg-q1_vect(1))-derivative_out);

    %system formulas
    system_states = veccat(plant_states, unit_delay_state, derivative_state, filter_state);
    system_noise_next = Function('system_noise_next',{system_states,system_noise_controls,params_sym}, {veccat( ...
        q1_vect, unit_delay_state_next, derivative_state_next, filter_state_next ...
        )},{'x','u','m'},{'x_next'}); %next state of system
    system_noise_vir_out = Function('system_noise_vir_out', {system_states,system_noise_controls,params_sym}, {sat_vir_out}); %system output: control voltage of plant


    %% [SYSTEM+newcontroller]
    input_noise = MX.sym('input_noise');
    input_qg = MX.sym('input_qg');

    %plant input    
    unit_delay_state = MX.sym('unit_delay_state');
    system_newcontroller_controls = veccat(input_qg, input_noise);
    %sat_vir_out = clamping_formula(unit_delay_state);
    sat_vir_out = unit_delay_state; %disabling clamping on vir

    %plant output
    q1_vect = plant_next(plant_states,sat_vir_out,params_sym)+input_noise;

    %filter block
    filter_state = MX.sym('filter_state',2,1);
    %y[k] = 2 u[k] - 2 u[k-1] + 0.998y[k-1]
    filter_out = 2*q1_vect(1)-2*filter_state(1)+0.998002030680999*filter_state(2);
    filter_state_next = [ q1_vect(1); filter_out ];    

    %unit delay
    unit_delay_state_next = kv * (kp*(input_qg-q1_vect(1))-filter_out);

    %system formulas
    system_newcontroller_states = veccat(plant_states, unit_delay_state, filter_state);
    system_newcontroller_next = Function('system_newcontroller_next', ...
        {system_newcontroller_states,system_newcontroller_controls,params_sym}, ...
        {veccat(q1_vect, unit_delay_state_next, filter_state_next )}, ...
        {'x','u','m'},{'x_next'}); %next state of system
    system_newcontroller_vir_out = Function('system_newcontroller_vir_out', {system_newcontroller_states,system_newcontroller_controls,params_sym}, {sat_vir_out}); %system output: control voltage of plant
