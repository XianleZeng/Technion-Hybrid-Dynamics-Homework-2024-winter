function [phi, phi_d, phi_dd]=angles_input(t)
% ANGLES_INPUT    Calculate time profiles of input joint angles and their derivations
% Input - time t 
% Output - vectors of phi(t), phi_d(t), phi_dd(t) 
% Xianle Zeng
% 10-Dec-2024 20:04:49

[omega,alfa,Psi,beta]=control_parameters;

if (t <= 2*pi/omega) 
	% phi 
	phi=zeros(2,1);
	phi(1)=-(omega^5*t^3*(alfa + beta*sin(Psi - omega*t))*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(32*pi^5);
	phi(2)=(omega^5*t^3*(alfa + beta*sin(Psi + omega*t))*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(32*pi^5);
	% phi_d 
	phi_d=zeros(2,1);
	phi_d(1)=(beta*omega^6*t^3*cos(Psi - omega*t)*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(32*pi^5) - (3*omega^5*t^2*(alfa + beta*sin(Psi - omega*t))*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(32*pi^5) - (omega^5*t^3*(12*t - (30*pi)/omega)*(alfa + beta*sin(Psi - omega*t)))/(32*pi^5);
	phi_d(2)=(omega^5*t^3*(12*t - (30*pi)/omega)*(alfa + beta*sin(Psi + omega*t)))/(32*pi^5) + (3*omega^5*t^2*(alfa + beta*sin(Psi + omega*t))*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(32*pi^5) + (beta*omega^6*t^3*cos(Psi + omega*t)*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(32*pi^5);
	% phi_dd 
	phi_dd=zeros(2,1);
	phi_dd(1)=(beta*omega^6*t^3*cos(Psi - omega*t)*(12*t - (30*pi)/omega))/(16*pi^5) - (3*omega^5*t^2*(12*t - (30*pi)/omega)*(alfa + beta*sin(Psi - omega*t)))/(16*pi^5) - (3*omega^5*t*(alfa + beta*sin(Psi - omega*t))*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(16*pi^5) - (3*omega^5*t^3*(alfa + beta*sin(Psi - omega*t)))/(8*pi^5) + (3*beta*omega^6*t^2*cos(Psi - omega*t)*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(16*pi^5) + (beta*omega^7*t^3*sin(Psi - omega*t)*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(32*pi^5);
	phi_dd(2)=(3*omega^5*t^3*(alfa + beta*sin(Psi + omega*t)))/(8*pi^5) + (3*omega^5*t^2*(12*t - (30*pi)/omega)*(alfa + beta*sin(Psi + omega*t)))/(16*pi^5) + (3*omega^5*t*(alfa + beta*sin(Psi + omega*t))*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(16*pi^5) + (beta*omega^6*t^3*cos(Psi + omega*t)*(12*t - (30*pi)/omega))/(16*pi^5) + (3*beta*omega^6*t^2*cos(Psi + omega*t)*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(16*pi^5) - (beta*omega^7*t^3*sin(Psi + omega*t)*(6*t^2 + (40*pi^2)/omega^2 - (30*t*pi)/omega))/(32*pi^5);
else 
	% phi 
	phi=zeros(2,1);
	phi(1)=- alfa - beta*sin(Psi - omega*t);
	phi(2)=alfa + beta*sin(Psi + omega*t);
	% phi_d 
	phi_d=zeros(2,1);
	phi_d(1)=beta*omega*cos(Psi - omega*t);
	phi_d(2)=beta*omega*cos(Psi + omega*t);
	% phi_dd 
	phi_dd=zeros(2,1);
	phi_dd(1)=beta*omega^2*sin(Psi - omega*t);
	phi_dd(2)=-beta*omega^2*sin(Psi + omega*t);
end 
