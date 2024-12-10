function theta_dd = theta_dd_H(q, q_d,t)
%THETA_DD Summary of this function goes here
%   Detailed explanation goes here

[m, l, g] = model_params;
x=q(1); y=q(2); th=q(3); phi1=q(4); phi2=q(5);
x_d=q_d(1); y_d=q_d(2); th_d=q_d(3); phi_d1=q_d(4); phi_d2=q_d(5);

%% Get the input angle (shape state)
[phi, phi_d, phi_dd]=angles_input(t);
phi_2 = phi(2);
phi_d2 = phi_d(2);
phi_dd2 = phi_dd(2);
phi_1 = phi(1);
phi_d1 = phi_d(1);
phi_dd1 = phi_dd(1);

%% Torque 
%
qb = q(1:3)';
qb_d = q_d(1:3)';
[qb_dd, tau]=dyn_sol(qb,qb_d,t); 
tau_1 = tau(1); tau_2 = tau(2); 
x_dd = qb_dd(1);
y_dd = qb_dd(2);

theta_dd = (tau_2 - m*g*l*cos(th + phi_2)  - (4*l^2*m/3)*phi_dd2 - l^2*m*sin(phi_2)*th_d^2 - l*m*cos(phi_2 + th)*y_dd + l*m*sin(phi_2 + th)*x_dd)/((4*l^2*m/3) + l^2*m*cos(phi_2));
end