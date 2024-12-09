function [p_cm, v_cm]=center_of_mass(q, q_d)
% CENTER_OF_MASS    center of mass of the cat
% Input - values of generalized coordinate and their velocity q, q_dot 
% Output - position of the center of mass p_cm, velocity of the center of mass v_cm 
% Xianle Zeng
% 09-Dec-2024 11:38:33

[m, l, g]=model_params;

x=q(1); y=q(2); th=q(3); phi1=q(4); phi2=q(5);
x_d=q_d(1); y_d=q_d(2); th_d=q_d(3); phi_d1=q_d(4); phi_d2=q_d(5);


% p_cm
p_cm=zeros(2,1);
p_cm(1)=(m*x + m*(x + l*(cos(phi2 + th) + cos(th))) + m*(x - l*(cos(phi1 + th) + cos(th))))/(3*m);
p_cm(2)=(m*y + m*(y + l*(sin(phi2 + th) + sin(th))) + m*(y - l*(sin(phi1 + th) + sin(th))))/(3*m);

% v_cm
v_cm=zeros(2,1);
v_cm(1)=x_d + (th_d*(l*m*(sin(phi1 + th) + sin(th)) - l*m*(sin(phi2 + th) + sin(th))))/(3*m) + (l*phi_d1*sin(phi1 + th))/3 - (l*phi_d2*sin(phi2 + th))/3;
v_cm(2)=y_d - (th_d*(l*m*(cos(phi1 + th) + cos(th)) - l*m*(cos(phi2 + th) + cos(th))))/(3*m) - (l*phi_d1*cos(phi1 + th))/3 + (l*phi_d2*cos(phi2 + th))/3;
end 
