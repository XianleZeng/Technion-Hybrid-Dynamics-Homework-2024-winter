function [q_dd, lambda]=dyn_sol_stick(q,q_d)
% DYN_SOL    Dynamic solution
% Xianle Zeng

%% Get the dynamic model H*[ddq; lambda] = K
[M,B,G,W,W_d]=dynamics_mat(q, q_d);

H = [M, -W.';
       W, zeros(2)];

K = [- B - G;
        -W_d*q_d];

dX = H\K;  

q_dd = dX(1:4);
lambda = dX(5:6);

end

