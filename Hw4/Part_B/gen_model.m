
clear

filename = strrep(mfilename,'gen_model','model_save');

if isempty(dir([filename,'.mat']))
    clear 

    filename = strrep(mfilename,'gen_model','model_save');

    %% variable declarations
    %
    syms x y theta phi real
    syms x_d y_d theta_d phi_d real
    syms x_dd y_dd theta_dd phi_dd real
    syms omega_0 v_0 theta_0 gamma beta real
    syms m_1 m_2 l h J_1 J_2 R real
    syms g real
    
    %% generalized coordinates
    %
    q = [x; y; theta; phi];
    dq = [x_d; y_d; theta_d; phi_d];
    ddq = [x_dd; y_dd; theta_dd; phi_dd];
    
    N=max(size(q));
    
    %% position, velocity of the center of mass
    %
    p_cm_horse = [x + h*sin(theta); y - h*cos(theta)];
    v_cm_horse = jacobian(p_cm_horse, q)*dq;
    p_cm_pendulum = [x + l*sin(phi+theta); y - l*cos(phi+theta)];
    v_cm_pendulum = jacobian(p_cm_pendulum, q)*dq;
    m = m_1 + m_2;
    p_cm_sys = (m_1*p_cm_horse + m_2*p_cm_pendulum)/m;
    v_cm_sys = jacobian(p_cm_sys, q)*dq;
    
    %% kinetic energy and ﻿potiential energy 
    %
    KE = (m_1/2)*(v_cm_horse'*v_cm_horse) + (m_2/2)*(v_cm_pendulum'*v_cm_pendulum) + (1/2)*J_1*theta_d^2 + (1/2)*J_2*(phi_d+theta_d)^2;
    PE = m*g*p_cm_sys(2);
    
    %% Constraints 
    %
    p_contact_point = [x; y - R];  % Contact point
    v_relative = cross([0; 0; theta_d],[0;-R;0]);
    v_contact_point = [x_d; y_d] + v_relative(1:2);
    h_ = p_contact_point(2);
    f_ = v_contact_point(1);
    
    W_h = jacobian(h_, q);
    W_f = jacobian(f_, dq);
    W = [W_f; W_h];
    
    num_constraints = 2;
    num_dof = 4;
    
    W_d = sym(zeros(num_constraints, num_dof));
    for i = 1:num_constraints
        for j = 1:num_dof
                W_d(i, j) = jacobian(W(i, j), q)*dq;
        end
    end
    
    %% M*ddq + C*dq + G = F + W^T*lambda
    % M
    
    M=simplify(jacobian(jacobian(KE,dq).',dq));
    
    %%
    syms C
    for k=1:N
        for j=1:N
            C(k,j)=0;
            for i=1:N
	            C(k,j)=C(k,j)+1/2*(diff(M(k,j),q(i)) + diff(M(k,i),q(j)) - diff(M(i,j),q(k)))*dq(i);
            end
        end
    end
    B = C*dq;
    
    %%
    %
    G = simplify(jacobian(PE, q))';
    
    %%
    
    lambda = simplify(inv(W*W.')*W*(M*ddq+B+G));


    %%
    %
	save(filename);
else
	clear
	load(strrep(mfilename,'gen_model','model_save'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Output functions for use in simulation
%%
%
%% First, output model to a file called
%%
%%    dynamics_mat.m
%%
%

%% Output file header
%
fcn_name='dynamics_mat';
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function [M,B,G,W,W_d]=%s',fcn_name);
fprintf(fid,'(q, q_d, damping)\n');
fprintf(fid,'%% %s    Model of three-link robot cat.\n',...
	upper(fcn_name));

fprintf(fid,'%% Xianle Zeng\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[m_1, m_2, l, h, J_1, J_2, R, g, mu]=model_params();\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'x=q(1); y=q(2); theta=q(3); phi=q(4); \n');
fprintf(fid,'x_d=q_d(1); y_d=q_d(2); theta_d=q_d(3); phi_d=q_d(4); \n\n');

%% Model output
%
fprintf(fid,'%% M matrix\n');
fprintf(fid,'M=zeros(%s);\n',num2str(N));
for k=1:N
	for j=1:N
		if M(k,j)~=0
			ttt=char(M(k,j));
			fprintf(fid,'M(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% B matrix\n');
fprintf(fid,'B=zeros(%s,1);\n',num2str(N));
for k=1:N
	if B(k)~=0
		ttt=char(B(k));
		fprintf(fid,'B(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'\n%% G matrix\n');
fprintf(fid,'G=zeros(%s,1);\n',num2str(N));
for k=1:N
	if G(k)~=0
		ttt=char(G(k));
		fprintf(fid,'G(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'\n%% W matrix\n');
fprintf(fid,'W=zeros(%s, %s);\n',num2str(num_constraints), num2str(N));
for k = 1:num_constraints
    for j = 1:N
		if W(k,j)~=0
			ttt=char(W(k,j));
			fprintf(fid,'W(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
    end
end


fprintf(fid,'\n%% W_d matrix\n');
fprintf(fid,'W_d=zeros(%s, %s);\n',num2str(num_constraints), num2str(N));
for k = 1:num_constraints
    for j = 1:N
		if W_d(k,j)~=0
			ttt=char(W_d(k,j));
			fprintf(fid,'W_d(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
    end
end

