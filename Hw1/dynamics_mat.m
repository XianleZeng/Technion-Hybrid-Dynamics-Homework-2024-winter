function [M,B,G]=dynamics_mat(q, q_d)
% DYNAMICS_MAT    Model of three-link robot cat.
% Input - values of generalized coordinate and their velocity q, q_dot 
% Output - matrices/vectors M(q), B(q, q_dot), G(q) for the dynamics equation of motion 
% Xianle Zeng
% 27-Nov-2024 19:50:49

[m, l, g]=model_params;

x=q(1); y=q(2); th=q(3); phi1=q(4); phi2=q(5);
x_d=q_d(1); y_d=q_d(2); th_d=q_d(3); phi_d1=q_d(4); phi_d2=q_d(5);

% M matrix
M=zeros(5);
M(1,1)=3*m;
M(1,3)=l*m*(sin(phi1 + th) - sin(phi2 + th));
M(1,4)=l*m*sin(phi1 + th);
M(1,5)=-l*m*sin(phi2 + th);
M(2,2)=3*m;
M(2,3)=-l*m*(cos(phi1 + th) - cos(phi2 + th));
M(2,4)=-l*m*cos(phi1 + th);
M(2,5)=l*m*cos(phi2 + th);
M(3,1)=l*m*(sin(phi1 + th) - sin(phi2 + th));
M(3,2)=-l*m*(cos(phi1 + th) - cos(phi2 + th));
M(3,3)=2*l^2*m*(cos(phi1) + cos(phi2) + 2);
M(3,4)=l^2*m*(cos(phi1) + 1);
M(3,5)=l^2*m*(cos(phi2) + 1);
M(4,1)=l*m*sin(phi1 + th);
M(4,2)=-l*m*cos(phi1 + th);
M(4,3)=l^2*m*(cos(phi1) + 1);
M(4,4)=l^2*m;
M(5,1)=-l*m*sin(phi2 + th);
M(5,2)=l*m*cos(phi2 + th);
M(5,3)=l^2*m*(cos(phi2) + 1);
M(5,5)=l^2*m;

% B matrix
B=zeros(5,1);
B(1)=th_d*(l*m*phi_d1*cos(phi1 + th) - l*m*phi_d2*cos(phi2 + th) + l*m*th_d*(cos(phi1 + th) - cos(phi2 + th))) + phi_d1*(l*m*phi_d1*cos(phi1 + th) + l*m*th_d*cos(phi1 + th)) - phi_d2*(l*m*phi_d2*cos(phi2 + th) + l*m*th_d*cos(phi2 + th));
B(2)=th_d*(l*m*phi_d1*sin(phi1 + th) - l*m*phi_d2*sin(phi2 + th) + l*m*th_d*(sin(phi1 + th) - sin(phi2 + th))) + phi_d1*(l*m*phi_d1*sin(phi1 + th) + l*m*th_d*sin(phi1 + th)) - phi_d2*(l*m*phi_d2*sin(phi2 + th) + l*m*th_d*sin(phi2 + th));
B(3)=- phi_d1*(l^2*m*phi_d1*sin(phi1) + l^2*m*th_d*sin(phi1)) - th_d*(l^2*m*phi_d1*sin(phi1) + l^2*m*phi_d2*sin(phi2)) - phi_d2*(l^2*m*phi_d2*sin(phi2) + l^2*m*th_d*sin(phi2));
B(4)=l^2*m*th_d^2*sin(phi1);
B(5)=l^2*m*th_d^2*sin(phi2);

% G matrix
G=zeros(5,1);
G(2)=-3*g*m;
G(3)=g*l*m*(cos(phi1 + th) + cos(th)) - g*l*m*(cos(phi2 + th) + cos(th));
G(4)=g*l*m*cos(phi1 + th);
G(5)=-g*l*m*cos(phi2 + th);
