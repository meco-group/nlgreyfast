function [dx, y] = emps_plant_atan_m(t, x, u, M1, Fv1, Fc1, OF1, varargin)
    % Function to define ODEs of the EMPS system for use with nlgreyest (System Identification Toolbox),
    % and for converting it to MEX.
    % Test call for MATLAB Coder: emps_plant_atan_m(0, [1;2], 2, 3, 4, 5, 6, {}) 
    % and specify a 0x0 cell inside an 1x1 cell for varargin...
    q=x(1); dq=x(2);
    gtau=35.150651882485469;
    %ddq = (1/M1)*(gtau*u - OF1 - sign(dq)*(Fv1*abs(dq)+Fc1)); %formula in the original paper
    ddq = (1/M1)*(gtau*u - OF1 - (Fc1*tanh(10000/Fc1*dq)+Fv1*dq)); %formula that we use, and it provides lower relative error on validation dataset with the parameters from the paper
    dx = [dq;ddq]; % state equations
    y = [q]; % output equations
end
