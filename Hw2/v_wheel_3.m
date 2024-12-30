function v_wheel_3=v_wheel_3(q, q_d)
% Xianle Zeng
% 30-Dec-2024 20:32:50

[m, l, I_c, d, w, b, c]=model_params;

x=q(1); y=q(2); theta=q(3); phi=q(4);
x_d=q_d(1); y_d=q_d(2); theta_d=q_d(3); phi_d=q_d(4);


% v_wheel_3
v_wheel_3=zeros(2,1);
v_wheel_3(1)=x_d - theta_d*(b*sin(phi + theta) + l*sin(theta)) - b*phi_d*sin(phi + theta);
v_wheel_3(2)=y_d + theta_d*(b*cos(phi + theta) + l*cos(theta)) + b*phi_d*cos(phi + theta);
