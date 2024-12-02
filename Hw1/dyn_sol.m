function [qb_dd, tau]=dyn_sol(qb,qb_d,t)
% DYN_SOL    Dynamic solution
% Input - time t, values of qb(t), qb_dot(t), 
% Output - body acceleration and joint torques tau(t), where tau = [tau1, tau2] 
% Xianle Zeng


%% Get the input angle (shape state)
[phi, phi_d, phi_dd]=angles_input(t);
qs = phi;
dqs = phi_d; 
ddqs = phi_dd;

q = [qb; qs];
q_d = [qb_d; dqs];

b = length(qb);
s = length(qs);

%% Get the dynamic model 
[M,B,G]=dynamics_mat(q, q_d);

Mbb = M(1:b, 1:b);
Mbs = M(1:b, b+1:end);
Mss = M(b+1:end, b+1:end);

Bb = B(1:b);
Bs = B(b+1:end);

Gb = G(1:b);
Gs = G(b+1:end);

A = [Mbb, zeros(b,s);
       Mbs', -eye(s,s)];

B = -[Mbs*ddqs + Bb + Gb;
        Mss*ddqs + Bs + Gs];

dX = A\B;

qb_dd = dX(1:3);
tau= dX(4:5);

end
