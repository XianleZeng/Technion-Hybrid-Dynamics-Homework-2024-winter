function [M,B,G,W,W_d]=dynamics_mat(q, q_d, damping)
% DYNAMICS_MAT    Model of three-link robot cat.
% Xianle Zeng
% 11-Feb-2025 23:20:02

[m_1, m_2, l, h, J_1, J_2, R, g, mu]=model_params();

x=q(1); y=q(2); theta=q(3); phi=q(4); 
x_d=q_d(1); y_d=q_d(2); theta_d=q_d(3); phi_d=q_d(4); 

% M matrix
M=zeros(4);
M(1,1)=m_1 + m_2;
M(1,3)=h*m_1*cos(theta) + l*m_2*cos(phi + theta);
M(1,4)=l*m_2*cos(phi + theta);
M(2,2)=m_1 + m_2;
M(2,3)=h*m_1*sin(theta) + l*m_2*sin(phi + theta);
M(2,4)=l*m_2*sin(phi + theta);
M(3,1)=h*m_1*cos(theta) + l*m_2*cos(phi + theta);
M(3,2)=h*m_1*sin(theta) + l*m_2*sin(phi + theta);
M(3,3)=J_1 + J_2 + h^2*m_1 + l^2*m_2;
M(3,4)=J_2 + l^2*m_2;
M(4,1)=l*m_2*cos(phi + theta);
M(4,2)=l*m_2*sin(phi + theta);
M(4,3)=J_2 + l^2*m_2;
M(4,4)=J_2 + l^2*m_2;

% B matrix
B=zeros(4,1);
B(1)=- theta_d*(theta_d*(h*m_1*sin(theta) + l*m_2*sin(phi + theta)) + l*m_2*phi_d*sin(phi + theta)) - phi_d*(l*m_2*phi_d*sin(phi + theta) + l*m_2*theta_d*sin(phi + theta));
B(2)=theta_d*(theta_d*(h*m_1*cos(theta) + l*m_2*cos(phi + theta)) + l*m_2*phi_d*cos(phi + theta)) + phi_d*(l*m_2*phi_d*cos(phi + theta) + l*m_2*theta_d*cos(phi + theta));

% G matrix
G=zeros(4,1);
G(2)=g*(m_1 + m_2);
G(3)=g*(h*m_1*sin(theta) + l*m_2*sin(phi + theta));
G(4)=g*l*m_2*sin(phi + theta);

% W matrix
W=zeros(2, 4);
W(1,1)=1;
W(1,3)=R;
W(2,2)=1;

% W_d matrix
W_d=zeros(2, 4);
