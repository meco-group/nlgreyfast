function [one_step, states, controls, params] = nlsim_one_step(ode, fs, N_steps_per_sample)

import casadi.*
dt = 1/fs/N_steps_per_sample;
states = ode.mx_in{1};
controls = ode.mx_in{2};
params = ode.mx_in{3};

% Build an integrator for this system: Runge Kutta 4 integrator
k1 = ode(states,controls,params);
k2 = ode(states+dt/2.0*k1,controls,params);
k3 = ode(states+dt/2.0*k2,controls,params);
k4 = ode(states+dt*k3,controls,params);
states_final = states+dt/6.0*(k1+2*k2+2*k3+k4);

% Create a function that simulates one step propagation in a sample
one_step = Function('one_step',{states, controls, params},{states_final});

