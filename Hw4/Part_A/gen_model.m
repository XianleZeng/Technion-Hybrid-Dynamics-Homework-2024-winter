clear

%% variable declarations
%
syms x y theta real
syms x_d y_d theta_d real
syms x_dd y_dd theta_dd real
syms omega_0 v_0 theta_0 gamma beta real
syms m l I_c real
syms g real

%% generalized coordinates
%
q = [x; y; theta];
dq = [x_d; y_d; theta_d];
ddq = [x_dd; y_dd; theta_dd];

N=max(size(q));

%% position, velocity of the center of mass
%
p_cm = [x; y];
v_cm = [x_d; y_d];

%% kinetic energy and ﻿potiential energy 
%
KE = simplify(m/2*(v_cm'*v_cm) + (1/2)*I_c*theta_d^2);
PE = m*g*y;

%% Holonomic Constraints 
%
h = y - l*cos(theta);
W = jacobian(h, q);

num_constraints = 1;
num_dof = 3;

W_d = sym(zeros(num_constraints, num_dof));
for i = 1:num_constraints
    for j = 1:num_dof
            W_d(i, j) = jacobian(W(i, j), q)*dq;
    end
end

%% M*ddq + C*dq + G = F + W^T*lambda
% M

M=simplify(jacobian(jacobian(KE,dq).',dq));

%%
syms C
for k=1:N
    for j=1:N
        C(k,j)=0;
        for i=1:N
	        C(k,j)=C(k,j)+1/2*(diff(M(k,j),q(i)) + diff(M(k,i),q(j)) - diff(M(i,j),q(k)))*dq(i);
        end
    end
end
B = C*dq;

%%
%
G = simplify(jacobian(PE, q))';

%%
% 
lambda_th_d = simplify(inv(W*inv(M)*W.')*(W*inv(M)*(B + G) - W_d*dq));
eq = subs(KE, {x_d, y_d}, {v_0, -l*sin(theta)*theta_d}) + subs(PE, {y}, {l*cos(theta)}) - subs(KE, {theta_d, x_d, y_d}, {omega_0, v_0, 0}) - m*g*l*cos(theta_0);
theta_d_sol = solve(eq, theta_d);
lambda_th = simplify(subs(lambda_th_d, theta_d, theta_d_sol(1)));
subs(lambda_th, {omega_0, I_c}, {gamma*g/l, beta*m*l^2})