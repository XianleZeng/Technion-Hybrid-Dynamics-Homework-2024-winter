clear;
close all;

%% Solve DAE
%
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
tspan = 0:0.05:100;  
X0 = [0; 0; 0; 0; 0; 0; 0; 0];             
[t, X] = ode45(@dyn_sol, tspan, X0, options); 
