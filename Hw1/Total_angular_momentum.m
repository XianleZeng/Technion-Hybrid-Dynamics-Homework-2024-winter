function H_total=Total_angular_momentum(q, q_d, t)
% TOTAL_ANGULAR_MOMENTUM    center of mass of the cat
% Input - values of generalized coordinate and their velocity q, q_dot 
% Output - Total Angular Momentum of the System H_tatol 
% Xianle Zeng
% 26-Dec-2024 11:27:13

[m, l, g]=model_params;

x=q(1); y=q(2); th=q(3); phi1=q(4); phi2=q(5);
x_d=q_d(1); y_d=q_d(2); th_d=q_d(3); phi_d1=q_d(4); phi_d2=q_d(5);

qb = [x; y; th]; qb_d = [x_d; y_d; th_d];
[qb_dd, tau]=dyn_sol(qb,qb_d,t); 
tau1 = tau(1); tau2 = tau(2); 

% H_total
H_total=m*(((m*y + m*(y + l*(sin(phi2 + th) + sin(th))) + m*(y - l*(sin(phi1 + th) + sin(th))))/(3*m) - y + l*(sin(phi1 + th) + sin(th)))*(l*th_d*(sin(phi1 + th) + sin(th)) - (th_d*(l*m*(sin(phi1 + th) + sin(th)) - l*m*(sin(phi2 + th) + sin(th))))/(3*m) + (2*l*phi_d1*sin(phi1 + th))/3 + (l*phi_d2*sin(phi2 + th))/3) + ((m*x + m*(x + l*(cos(phi2 + th) + cos(th))) + m*(x - l*(cos(phi1 + th) + cos(th))))/(3*m) - x + l*(cos(phi1 + th) + cos(th)))*(l*th_d*(cos(phi1 + th) + cos(th)) - (th_d*(l*m*(cos(phi1 + th) + cos(th)) - l*m*(cos(phi2 + th) + cos(th))))/(3*m) + (2*l*phi_d1*cos(phi1 + th))/3 + (l*phi_d2*cos(phi2 + th))/3)) + m*((y - (m*y + m*(y + l*(sin(phi2 + th) + sin(th))) + m*(y - l*(sin(phi1 + th) + sin(th))))/(3*m) + l*(sin(phi2 + th) + sin(th)))*(l*th_d*(sin(phi2 + th) + sin(th)) + (th_d*(l*m*(sin(phi1 + th) + sin(th)) - l*m*(sin(phi2 + th) + sin(th))))/(3*m) + (l*phi_d1*sin(phi1 + th))/3 + (2*l*phi_d2*sin(phi2 + th))/3) + (x - (m*x + m*(x + l*(cos(phi2 + th) + cos(th))) + m*(x - l*(cos(phi1 + th) + cos(th))))/(3*m) + l*(cos(phi2 + th) + cos(th)))*(l*th_d*(cos(phi2 + th) + cos(th)) + (th_d*(l*m*(cos(phi1 + th) + cos(th)) - l*m*(cos(phi2 + th) + cos(th))))/(3*m) + (l*phi_d1*cos(phi1 + th))/3 + (2*l*phi_d2*cos(phi2 + th))/3)) + m*((y - (m*y + m*(y + l*(sin(phi2 + th) + sin(th))) + m*(y - l*(sin(phi1 + th) + sin(th))))/(3*m))*((th_d*(l*m*(sin(phi1 + th) + sin(th)) - l*m*(sin(phi2 + th) + sin(th))))/(3*m) + (l*phi_d1*sin(phi1 + th))/3 - (l*phi_d2*sin(phi2 + th))/3) + (x - (m*x + m*(x + l*(cos(phi2 + th) + cos(th))) + m*(x - l*(cos(phi1 + th) + cos(th))))/(3*m))*((th_d*(l*m*(cos(phi1 + th) + cos(th)) - l*m*(cos(phi2 + th) + cos(th))))/(3*m) + (l*phi_d1*cos(phi1 + th))/3 - (l*phi_d2*cos(phi2 + th))/3)) + (l^2*m*th_d)/3 + (l^2*m*(phi_d1 + th_d))/3 + (l^2*m*(phi_d2 + th_d))/3;
end 
