function [qp_dd, lambda]=dyn_sol(qp,dqp,t)
% DYN_SOL    Dynamic solution
% Xianle Zeng


%% Get the input angle (shape state)
[phi, phi_d, phi_dd]=angles_input(t);
qa = phi;
dqa = phi_d; 
ddqa = phi_dd;

q = [qp; qa];
dq = [dqp; dqa];

n_p = length(qp);
n_a = length(qa);

%% Get the dynamic model 
[M,B,W,W_d]=dynamics_mat(q, dq);

Mpp = M(1:n_p, 1:n_p);
Mpa = M(1:n_a, n_a+1:end);
Maa = M(n_a+1:end, n_a+1:end);

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

