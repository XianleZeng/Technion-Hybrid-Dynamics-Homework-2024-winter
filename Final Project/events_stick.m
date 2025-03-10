function [value, isterminal, direction] = events_stick(t,X)

[m, m_h, l, I_c, g, alpha, mu]=model_params();

q = X(1:4);  q_d = X(5:8);
x=q(1); y=q(2); theta_1=q(3); theta_2=q(4); 
dx=q_d(1); dy=q_d(2); dtheta_1=q_d(3); dtheta_2=q_d(4); 

[~, lambda] = dyn_sol_stick(q,q_d);  % A*[q_dd; lambda] = B

lambda_t = lambda(1);
lambda_n = lambda(2);

val_minus = lambda_t - mu*lambda_n;
ister_minus = 1;
dir_minus = 1;

val_plus = lambda_t + mu*lambda_n;
ister_plus = 1;
dir_plus = -1;

H = y + 2*l*cos(theta_1) - 2*l*cos(theta_2);
H_d = dy - 2*l*sin(theta_1)*dtheta_1 +2*l*sin(theta_2)*dtheta_2;


val_failure = H;
ister_failure = 0;
dir_failure = -1;

if H_d < -0.1
    ister_failure = 1;
end

value = [val_plus; val_minus; val_failure];
isterminal =[ister_plus; ister_minus; ister_failure];
direction = [dir_plus; dir_minus; dir_failure];

end

