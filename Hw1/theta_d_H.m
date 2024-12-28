function th_d_H=theta_d_H(q, q_d)
% THETA_D_H    center of mass of the cat
% Xianle Zeng
% 28-Dec-2024 19:40:30

[m, l, g]=model_params;

x=q(1); y=q(2); th=q(3); phi1=q(4); phi2=q(5);
x_d=q_d(1); y_d=q_d(2); th_d=q_d(3); phi_d1=q_d(4); phi_d2=q_d(5);


% th_d_H
th_d_H=-(m*(((l*phi_d1*sin(phi1 + th))/3 - (l*phi_d2*sin(phi2 + th))/3)*(y - (m*y + m*(y + l*(sin(phi2 + th) + sin(th))) + m*(y - l*(sin(phi1 + th) + sin(th))))/(3*m)) + ((l*phi_d1*cos(phi1 + th))/3 - (l*phi_d2*cos(phi2 + th))/3)*(x - (m*x + m*(x + l*(cos(phi2 + th) + cos(th))) + m*(x - l*(cos(phi1 + th) + cos(th))))/(3*m))) + m*(((2*l*phi_d1*sin(phi1 + th))/3 + (l*phi_d2*sin(phi2 + th))/3)*((m*y + m*(y + l*(sin(phi2 + th) + sin(th))) + m*(y - l*(sin(phi1 + th) + sin(th))))/(3*m) - y + l*(sin(phi1 + th) + sin(th))) + ((2*l*phi_d1*cos(phi1 + th))/3 + (l*phi_d2*cos(phi2 + th))/3)*((m*x + m*(x + l*(cos(phi2 + th) + cos(th))) + m*(x - l*(cos(phi1 + th) + cos(th))))/(3*m) - x + l*(cos(phi1 + th) + cos(th)))) + m*(((l*phi_d1*cos(phi1 + th))/3 + (2*l*phi_d2*cos(phi2 + th))/3)*(x - (m*x + m*(x + l*(cos(phi2 + th) + cos(th))) + m*(x - l*(cos(phi1 + th) + cos(th))))/(3*m) + l*(cos(phi2 + th) + cos(th))) + ((l*phi_d1*sin(phi1 + th))/3 + (2*l*phi_d2*sin(phi2 + th))/3)*(y - (m*y + m*(y + l*(sin(phi2 + th) + sin(th))) + m*(y - l*(sin(phi1 + th) + sin(th))))/(3*m) + l*(sin(phi2 + th) + sin(th)))) + (l^2*m*phi_d1)/3 + (l^2*m*phi_d2)/3)/(m*((l*(cos(phi2 + th) + cos(th)) + (l*m*(cos(phi1 + th) + cos(th)) - l*m*(cos(phi2 + th) + cos(th)))/(3*m))*(x - (m*x + m*(x + l*(cos(phi2 + th) + cos(th))) + m*(x - l*(cos(phi1 + th) + cos(th))))/(3*m) + l*(cos(phi2 + th) + cos(th))) + (l*(sin(phi2 + th) + sin(th)) + (l*m*(sin(phi1 + th) + sin(th)) - l*m*(sin(phi2 + th) + sin(th)))/(3*m))*(y - (m*y + m*(y + l*(sin(phi2 + th) + sin(th))) + m*(y - l*(sin(phi1 + th) + sin(th))))/(3*m) + l*(sin(phi2 + th) + sin(th)))) + l^2*m + m*(((y - (m*y + m*(y + l*(sin(phi2 + th) + sin(th))) + m*(y - l*(sin(phi1 + th) + sin(th))))/(3*m))*(l*m*(sin(phi1 + th) + sin(th)) - l*m*(sin(phi2 + th) + sin(th))))/(3*m) + ((x - (m*x + m*(x + l*(cos(phi2 + th) + cos(th))) + m*(x - l*(cos(phi1 + th) + cos(th))))/(3*m))*(l*m*(cos(phi1 + th) + cos(th)) - l*m*(cos(phi2 + th) + cos(th))))/(3*m)) + m*((l*(sin(phi1 + th) + sin(th)) - (l*m*(sin(phi1 + th) + sin(th)) - l*m*(sin(phi2 + th) + sin(th)))/(3*m))*((m*y + m*(y + l*(sin(phi2 + th) + sin(th))) + m*(y - l*(sin(phi1 + th) + sin(th))))/(3*m) - y + l*(sin(phi1 + th) + sin(th))) + (l*(cos(phi1 + th) + cos(th)) - (l*m*(cos(phi1 + th) + cos(th)) - l*m*(cos(phi2 + th) + cos(th)))/(3*m))*((m*x + m*(x + l*(cos(phi2 + th) + cos(th))) + m*(x - l*(cos(phi1 + th) + cos(th))))/(3*m) - x + l*(cos(phi1 + th) + cos(th)))));
end 
