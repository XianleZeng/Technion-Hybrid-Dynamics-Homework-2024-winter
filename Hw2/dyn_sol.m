function [q_dd, lambda]=dyn_sol(t,X)
% DYN_SOL    Dynamic solution
% Xianle Zeng


%% Get tactuated torque
%
F_qa = input_torque(t);

%% read in params
%
m = 2;
x=X(1); y=X(2); theta=X(3); phi=X(4); 
x_d=X(5); y_d=X(6); theta_d=X(7); phi_d=X(8); 

q_p = [x; y; theta];
q_a = phi;
dq_p = [x_d; y_d; theta_d];
dq_a = phi_d;
q = [q_p; q_a];
dq = [dq_p; dq_a];

%% Get the dynamic model 
[M,Mpp,Mpa,Maa,B,Bp,Ba,W,W_d,Wp,Wa,Wp_d,Wa_d]=dynamics_mat(q, dq);

% H = [Mpp, Mpa, -Wp';
%        Mpa', Maa, -Wa';
%        Wp, Wa, zeros(m)];
H = [M, -W';
       W, zeros(2)];

% K = [-Bp - zeros(size(Bp));
%         F_qa - Ba - zeros(size(Ba));
%         -Wp_d*dq_p - Wa_d*dq_a];

F_q = [0; 0; 0; F_qa];
K = [F_q - B - zeros(size(B));
        -W_d*dq];

dX = H\K;

q_dd = zeros(8,1);
q_dd(1:4) = dq;
q_dd(5:8) = dX(1:4);
lambda= dX(5:6);

end
