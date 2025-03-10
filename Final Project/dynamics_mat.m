function [M,B,G,W,W_d]=dynamics_mat(q, q_d)
% DYNAMICS_MAT    Model of three-link robot cat.
% Xianle Zeng
% 11-Mar-2025 00:02:09

[m, m_h, l, I_c, g, alpha, mu]=model_params();

x=q(1); y=q(2); theta_1=q(3); theta_2=q(4); 
dx=q_d(1); dy=q_d(2); dtheta_1=q_d(3); dtheta_2=q_d(4); 

% M matrix
M=zeros(4);
M(1,1)=2*m + m_h;
M(1,3)=l*cos(theta_1)*(3*m + 2*m_h);
M(1,4)=l*m*cos(theta_2);
M(2,2)=2*m + m_h;
M(2,3)=-l*sin(theta_1)*(3*m + 2*m_h);
M(2,4)=l*m*sin(theta_2);
M(3,1)=l*cos(theta_1)*(3*m + 2*m_h);
M(3,2)=-l*sin(theta_1)*(3*m + 2*m_h);
M(3,3)=I_c + 5*l^2*m + 4*l^2*m_h;
M(3,4)=2*l^2*m*cos(theta_1 + theta_2);
M(4,1)=l*m*cos(theta_2);
M(4,2)=l*m*sin(theta_2);
M(4,3)=2*l^2*m*cos(theta_1 + theta_2);
M(4,4)=I_c + l^2*m;

% B matrix
B=zeros(4,1);
B(1)=- dtheta_1^2*l*sin(theta_1)*(3*m + 2*m_h) - dtheta_2^2*l*m*sin(theta_2);
B(2)=dtheta_2^2*l*m*cos(theta_2) - dtheta_1^2*l*cos(theta_1)*(3*m + 2*m_h);
B(3)=-2*dtheta_2^2*l^2*m*sin(theta_1 + theta_2);
B(4)=-2*dtheta_1^2*l^2*m*sin(theta_1 + theta_2);

% G matrix
G=zeros(4,1);
G(1)=-g*sin(alpha)*(2*m + m_h);
G(2)=g*cos(alpha)*(2*m + m_h);
G(3)=-g*l*sin(alpha + theta_1)*(3*m + 2*m_h);
G(4)=-g*l*m*sin(alpha - theta_2);

% W matrix
W=zeros(2, 4);
W(1,1)=1;
W(2,2)=1;

% W_d matrix
W_d=zeros(2, 4);
