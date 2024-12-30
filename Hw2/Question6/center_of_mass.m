function [r_cm, v_cm, a_cm]=center_of_mass(q, q_d, q_dd)
% CENTER_OF_MASS    center of mass of the cat
% Input - values of generalized coordinate and their velocity q, q_dot 
% Output - position of the center of mass p_cm, velocity of the center of mass v_cm 
% Xianle Zeng
% 30-Dec-2024 21:26:26

[m, l, I_c, d, w, b, c]=model_params;

x=q(1); y=q(2); theta=q(3); phi=q(4);
x_d=q_d(1); y_d=q_d(2); theta_d=q_d(3); phi_d=q_d(4);

x_dd=q_dd(1); y_dd=q_dd(2); theta_dd=q_dd(3); phi_dd=q_dd(4);


% r_cm
r_cm=zeros(2,1);
r_cm(1)=x + d*cos(theta);
r_cm(2)=y + d*sin(theta);

% v_cm
v_cm=zeros(2,1);
v_cm(1)=x_d - d*theta_d*sin(theta);
v_cm(2)=y_d + d*theta_d*cos(theta);

% a_cm
a_cm=zeros(2,1);
a_cm(1)=x_dd - d*theta_dd*sin(theta) - d*theta_d^2*cos(theta);
a_cm(2)=y_dd + d*theta_dd*cos(theta) - d*theta_d^2*sin(theta);
end 
