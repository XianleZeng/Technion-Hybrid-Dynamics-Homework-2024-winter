function [q_dd, lambda]=dyn_sol_Q3(t,q,q_d)
% DYN_SOL    Dynamic solution
% Xianle Zeng


%% Get actuated torque
%
F_qa = input_torque(t);


%% Get the dynamic model H*[ddq; lambda] = K
[M,Mpp,Mpa,Maa,B,Bp,Ba,G,Gp,Ga,W,W_d,Wp,Wa,Wp_d,Wa_d]=dynamics_mat(q, q_d, false);

H = [M, -W';
       W, zeros(2)];

F_q = [0; 0; 0; F_qa];
K = [F_q - B - G;
        -W_d*q_d];

dX = H\K;  

q_dd = dX(1:4);
lambda = dX(5:6);

end

