function [one_sample, all_samples] = nlsim(ode, fs, N, N_steps_per_sample)
% Discretize *ode* at *1/fs* time step, return both function *one_step* of the next state (of the discretized system) or *all_samples* of multiple next states (starting from a given initial state). Both shall be computationally efficient (see *expand()* and *mapaccum()* of CasADi).

import casadi.*
[one_step, states, controls, params] = nlsim_one_step(ode, fs, N_steps_per_sample);

X = states;
for i=1:N_steps_per_sample
    X = one_step(X, controls, params);
end

% Create a function that simulates all step propagation on a sample
one_sample = Function('one_sample',{states, controls, params}, {X}, {'x','u','p'},{'x_next'});

% speedup trick: expand into scalar operations
one_sample = one_sample.expand();

%% Simulating the system
if nargout > 1
    all_samples = one_sample.mapaccum('all_samples', N); %mapaccum is a bit faster than the for loop
end
