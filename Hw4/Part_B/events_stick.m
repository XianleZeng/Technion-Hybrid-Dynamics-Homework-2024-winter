function [value, isterminal, direction] = events_stick(t,X)

[m_1, m_2, l, h, J_1, J_2, R, g, mu]=model_params();

q = X(1:4);  q_d = X(5:8);
x=q(1); y=q(2); theta=q(3); phi=q(4); 
x_d=q_d(1); y_d=q_d(2); theta_d=q_d(3); phi_d=q_d(4); 

[~, lambda] = dyn_sol_stick(t,q,q_d);  % A*[q_dd; lambda] = B

lambda_t = lambda(1);
lambda_n = lambda(2);

val_minus = lambda_t - mu*lambda_n;
ister_minus = 1;
dir_minus = 1;

val_plus = lambda_t + mu*lambda_n;
ister_plus = 1;
dir_plus = -1;

value = [val_plus; val_minus];
isterminal =[ister_plus; ister_minus];
direction = [dir_plus; dir_minus];

end

