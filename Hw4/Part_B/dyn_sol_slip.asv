function [q_dd, lambda]=dyn_sol_slip(t,q,q_d,slip_dir)
% DYN_SOL    Dynamic solution
% Xianle Zeng

%% Get the dynamic model H*[ddq; lambda] = K
[M,B,G,W,W_d]=dynamics_mat(q, q_d);

[m_1, m_2, l, h, J_1, J_2, R, g, mu]=model_params();

W_t = W(1, :);
W_n = W(2, :);
Wd_t = W_d(1, :);
Wd_n = W_d(1, :);

if (slip_dir == 1)
    sigma = 1;
elseif (slip_dir == 2)
    sigma = -1;
end

H = [M, -(W_n - sigma*mu*W_t).';
       W_n, zeros(1)];

K = [- B - G;
        -Wd_n*q_d];

dX = H\K;  

q_dd = dX(1:4);
lambda = dX(5);

end
