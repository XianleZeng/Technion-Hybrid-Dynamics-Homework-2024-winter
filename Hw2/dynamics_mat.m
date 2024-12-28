function [M,Mpp,Mpa,Maa,B,Bp,Ba,W,W_d,Wp,Wa,Wp_d,Wa_d]=dynamics_mat(q, q_d)
% DYNAMICS_MAT    Model of three-link robot cat.
% Xianle Zeng
% 28-Dec-2024 21:44:28

[m, l, I_c, d, w, b, c]=model_params;

x=q(1); y=q(2); theta=q(3); phi=q(4); 
x_d=q_d(1); y_d=q_d(2); theta_d=q_d(3); phi_d=q_d(4); 

% M matrix
M=zeros(4);
M(1,1)=m;
M(1,3)=-d*m*sin(theta);
M(2,2)=m;
M(2,3)=d*m*cos(theta);
M(3,1)=-d*m*sin(theta);
M(3,2)=d*m*cos(theta);
M(3,3)=I_c + d^2*m;
% Mpp matrix
Mpp=zeros(3);
Mpp(1,1)=m;
Mpp(1,3)=-d*m*sin(theta);
Mpp(2,2)=m;
Mpp(2,3)=d*m*cos(theta);
Mpp(3,1)=-d*m*sin(theta);
Mpp(3,2)=d*m*cos(theta);
Mpp(3,3)=I_c + d^2*m;
% Maa matrix
Maa=zeros(1);
% Mpa matrix
Mpa=zeros(3, 1);

% B matrix
B=zeros(4,1);
B(1)=- c*(l*theta_d*sin(theta) - 3*x_d + b*phi_d*sin(phi + theta) + b*theta_d*sin(phi + theta)) - d*m*theta_d^2*cos(theta);
B(2)=c*(3*y_d + l*theta_d*cos(theta) + b*phi_d*cos(phi + theta) + b*theta_d*cos(phi + theta)) - d*m*theta_d^2*sin(theta);
B(3)=(c*(2*b^2*phi_d + 2*b^2*theta_d + 2*l^2*theta_d + theta_d*w^2 + 2*l*y_d*cos(theta) - 2*l*x_d*sin(theta) + 2*b*y_d*cos(phi + theta) - 2*b*x_d*sin(phi + theta) + 2*b*l*phi_d*cos(phi) + 4*b*l*theta_d*cos(phi)))/2;
B(4)=b*c*(b*phi_d + b*theta_d + y_d*cos(phi + theta) - x_d*sin(phi + theta) + l*theta_d*cos(phi));

% Bp matrix
Bp=zeros(3,1);
Bp(1)=- c*(l*theta_d*sin(theta) - 3*x_d + b*phi_d*sin(phi + theta) + b*theta_d*sin(phi + theta)) - d*m*theta_d^2*cos(theta);
Bp(2)=c*(3*y_d + l*theta_d*cos(theta) + b*phi_d*cos(phi + theta) + b*theta_d*cos(phi + theta)) - d*m*theta_d^2*sin(theta);
Bp(3)=(c*(2*b^2*phi_d + 2*b^2*theta_d + 2*l^2*theta_d + theta_d*w^2 + 2*l*y_d*cos(theta) - 2*l*x_d*sin(theta) + 2*b*y_d*cos(phi + theta) - 2*b*x_d*sin(phi + theta) + 2*b*l*phi_d*cos(phi) + 4*b*l*theta_d*cos(phi)))/2;

% Ba matrix
Ba=zeros(1,1);
Ba(1)=b*c*(b*phi_d + b*theta_d + y_d*cos(phi + theta) - x_d*sin(phi + theta) + l*theta_d*cos(phi));

% W matrix
W=zeros(2, 4);
W(1,1)=-sin(theta);
W(1,2)=cos(theta);
W(2,1)=-sin(phi + theta);
W(2,2)=cos(phi + theta);
W(2,3)=b + l*cos(phi);
W(2,4)=b;

% Wp matrix
Wp=zeros(2, 3);
Wp(1,1)=-sin(theta);
Wp(1,2)=cos(theta);
Wp(2,1)=-sin(phi + theta);
Wp(2,2)=cos(phi + theta);
Wp(2,3)=b + l*cos(phi);

% Wa matrix
Wa=zeros(2,1);
Wa(2)=b;

% W_d matrix
W_d=zeros(2, 4);
W_d(1,1)=-theta_d*cos(theta);
W_d(1,2)=-theta_d*sin(theta);
W_d(2,1)=- phi_d*cos(phi + theta) - theta_d*cos(phi + theta);
W_d(2,2)=- phi_d*sin(phi + theta) - theta_d*sin(phi + theta);
W_d(2,3)=-l*phi_d*sin(phi);

% Wp_d matrix
Wp_d=zeros(2, 3);
Wp_d(1,1)=-theta_d*cos(theta);
Wp_d(1,2)=-theta_d*sin(theta);
Wp_d(2,1)=- phi_d*cos(phi + theta) - theta_d*cos(phi + theta);
Wp_d(2,2)=- phi_d*sin(phi + theta) - theta_d*sin(phi + theta);
Wp_d(2,3)=-l*phi_d*sin(phi);

% Wa_d matrix
Wa_d=zeros(2,1);
