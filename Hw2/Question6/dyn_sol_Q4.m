function [qp_dd, F_qa, lambda]=dyn_sol_Q4(t,qp,qp_d,damping)
% DYN_SOL_Q4    Dynamic solution
% Xianle Zeng


%% Get input angle
%
[phi, phi_d, phi_dd] = input_angle(t);

%% read in params
%
m = 2;
x=qp(1); y=qp(2); theta=qp(3);
x_d=qp_d(1); y_d=qp_d(2); theta_d=qp_d(3); 

q_p = [x; y; theta];
q_a = phi;
dq_p = [x_d; y_d; theta_d];
dq_a = phi_d;
ddq_a = phi_dd;
q = [q_p; q_a];
dq = [dq_p; dq_a];

n_a = length(q_a);
n_p = length(q_p);


%% Get the dynamic model H*[ddq; lambda] = K
[M,Mpp,Mpa,Maa,B,Bp,Ba,G,Gp,Ga,W,W_d,Wp,Wa,Wp_d,Wa_d]=dynamics_mat(q, dq, damping);

H = [Mpp, Mpa, -Wp';
       Mpa', -eye(n_a, n_a), -Wa';
       Wp, zeros(m, n_a), zeros(m)];

K = -[Mpa*ddq_a + Bp + Gp;
       Maa*ddq_a + Ba + Ga;
       Wa*ddq_a + Wp_d*dq_p + Wa_d*dq_a];

dX = H\K;

qp_dd = dX(1:3);
F_qa = dX(4);
lambda= dX(5:6);

end

