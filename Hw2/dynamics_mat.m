function [M,Mpp,Mpa,Maa,B,Bp,Ba,W,W_d,Wp,Wa,Wp_d,Wa_d]=dynamics_mat(q, q_d)
% DYNAMICS_MAT    Model of three-link robot cat.
% Xianle Zeng
% 28-Dec-2024 20:42:00

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
B(1)=- c*(3*d*theta_d*sin(theta) - 3*x_d + b*phi_d*sin(phi + theta) + b*theta_d*sin(phi + theta)) - d*m*theta_d^2*cos(theta);
B(2)=c*(y_d + (phi_d + theta_d)*(l - d + b*cos(phi + theta)) + d*theta_d*cos(theta)) + 2*c*(y_d - d*theta_d + d*theta_d*cos(theta)) - d*m*theta_d^2*sin(theta);
B(3)=(c*(2*(w/2 - d*sin(theta))*(x_d + (theta_d*w)/2 - d*theta_d*sin(theta)) - 2*(d - d*cos(theta))*(y_d - d*theta_d + d*theta_d*cos(theta))))/2 + (c*(2*(y_d + (phi_d + theta_d)*(l - d + b*cos(phi + theta)) + d*theta_d*cos(theta))*(l - d + b*cos(phi + theta) + d*cos(theta)) + 2*(b*sin(phi + theta) + d*sin(theta))*(d*theta_d*sin(theta) - x_d + b*sin(phi + theta)*(phi_d + theta_d))))/2 - (c*(2*(d - d*cos(theta))*(y_d - d*theta_d + d*theta_d*cos(theta)) - 2*(w/2 + d*sin(theta))*((theta_d*w)/2 - x_d + d*theta_d*sin(theta))))/2;
B(4)=(c*(2*(y_d + (phi_d + theta_d)*(l - d + b*cos(phi + theta)) + d*theta_d*cos(theta))*(l - d + b*cos(phi + theta)) + 2*b*sin(phi + theta)*(d*theta_d*sin(theta) - x_d + b*sin(phi + theta)*(phi_d + theta_d))))/2;

% Bp matrix
Bp=zeros(3,1);
Bp(1)=- c*(3*d*theta_d*sin(theta) - 3*x_d + b*phi_d*sin(phi + theta) + b*theta_d*sin(phi + theta)) - d*m*theta_d^2*cos(theta);
Bp(2)=c*(y_d + (phi_d + theta_d)*(l - d + b*cos(phi + theta)) + d*theta_d*cos(theta)) + 2*c*(y_d - d*theta_d + d*theta_d*cos(theta)) - d*m*theta_d^2*sin(theta);
Bp(3)=(c*(2*(w/2 - d*sin(theta))*(x_d + (theta_d*w)/2 - d*theta_d*sin(theta)) - 2*(d - d*cos(theta))*(y_d - d*theta_d + d*theta_d*cos(theta))))/2 + (c*(2*(y_d + (phi_d + theta_d)*(l - d + b*cos(phi + theta)) + d*theta_d*cos(theta))*(l - d + b*cos(phi + theta) + d*cos(theta)) + 2*(b*sin(phi + theta) + d*sin(theta))*(d*theta_d*sin(theta) - x_d + b*sin(phi + theta)*(phi_d + theta_d))))/2 - (c*(2*(d - d*cos(theta))*(y_d - d*theta_d + d*theta_d*cos(theta)) - 2*(w/2 + d*sin(theta))*((theta_d*w)/2 - x_d + d*theta_d*sin(theta))))/2;

% Ba matrix
Ba=zeros(1,1);
Ba(1)=(c*(2*(y_d + (phi_d + theta_d)*(l - d + b*cos(phi + theta)) + d*theta_d*cos(theta))*(l - d + b*cos(phi + theta)) + 2*b*sin(phi + theta)*(d*theta_d*sin(theta) - x_d + b*sin(phi + theta)*(phi_d + theta_d))))/2;

% W matrix
W=zeros(2, 4);
W(1,1)=-sin(theta);
W(1,2)=cos(theta);
W(2,1)=-1;
W(2,3)=b*sin(phi + theta) + d*sin(theta);
W(2,4)=b*sin(phi + theta);

% Wp matrix
Wp=zeros(2, 3);
Wp(1,1)=-sin(theta);
Wp(1,2)=cos(theta);
Wp(2,1)=-1;
Wp(2,3)=b*sin(phi + theta) + d*sin(theta);

% Wa matrix
Wa=zeros(2,1);
Wa(2)=b*sin(phi + theta);

% W_d matrix
W_d=zeros(2, 4);
W_d(1,1)=-theta_d*cos(theta);
W_d(1,2)=-theta_d*sin(theta);
W_d(2,3)=theta_d*(b*cos(phi + theta) + d*cos(theta)) + b*phi_d*cos(phi + theta);
W_d(2,4)=b*phi_d*cos(phi + theta) + b*theta_d*cos(phi + theta);

% Wp_d matrix
Wp_d=zeros(2, 3);
Wp_d(1,1)=-theta_d*cos(theta);
Wp_d(1,2)=-theta_d*sin(theta);
Wp_d(2,3)=theta_d*(b*cos(phi + theta) + d*cos(theta)) + b*phi_d*cos(phi + theta);

% Wa_d matrix
Wa_d=zeros(2,1);
Wa_d(2)=b*phi_d*cos(phi + theta) + b*theta_d*cos(phi + theta);
