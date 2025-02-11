function [value, isterminal, direction] = events_stick(t,X,slip_dir)

persistent first_call
if isempty(first_call)
    first_call = true;
end

[m_1, m_2, l, h, J_1, J_2, R, g, mu]=model_params();

q = X(1:4);  q_d = X(5:8);
x=q(1); y=q(2); theta=q(3); phi=q(4); 
x_d=q_d(1); y_d=q_d(2); theta_d=q_d(3); phi_d=q_d(4); 

[~, lambda_n] = dyn_sol_slip(t,q,q_d,slip_dir);  % A*[q_dd; lambda] = B

if (slip_dir == 1)
    sigma = 1;
else
    sigma = -1;
end

vt = x_d - R*phi_d;

value = [vt; lambda_n];

if first_call
    value(1) = 1;  % 避免事件在 t0 直接触发
    first_call = false;
end

isterminal =[1; 1];
direction = [-sigma; -1];

end

