
clear

filename = strrep(mfilename,'gen_model','model_save');

if isempty(dir([filename,'.mat']))
    clear 

    filename = strrep(mfilename,'gen_model','model_save');

    syms x y theta_1 theta_2 real
    syms dx dy dtheta_1 dtheta_2 real
    syms ddx ddy ddtheta_1 ddtheta_2 real
    syms alpha m_h m g l I_c real
    
    %% generalized coordinates
    %
    q = [x; y; theta_1; theta_2];
    dq = [dx; dy; dtheta_1; dtheta_2];
    ddq = [ddx; ddy; ddtheta_1; ddtheta_2];
    N=max(size(q));
    
    %% position of masses in system (center of mass for each solid body)
    %
    r_stance = [x + l*sin(theta_1); y + l*cos(theta_1)];
    r_swing = [x + 2*l*sin(theta_1) + l*sin(theta_2); y + 2*l*cos(theta_1) - l*cos(theta_2)];
    r_hip = [x + 2*l*sin(theta_1); y + 2*l*cos(theta_1)];
    
    %% velocities of masses in system (center of mass for each solid body)
    %
    v_stance = jacobian(r_stance,q)*dq;
    v_swing = jacobian(r_swing,q)*dq;
    v_hip = jacobian(r_hip,q)*dq;
    
    %% kinetic energy of masses in system
    %
    KE_stance = simplify(m/2*(v_stance'*v_stance)) + (1/2)*I_c*dtheta_1^2;
    KE_swing = simplify(m/2*(v_swing'*v_swing)) + (1/2)*I_c*dtheta_2^2;
    KE_hip = simplify(m_h/2*(v_hip'*v_hip));
    KE = KE_stance+KE_swing+KE_hip;
    
    %% potential energy of masses in system
    %
    g_vec = [sin(alpha); -cos(alpha)]*g;
    PE_stance = r_stance.'*g_vec*m;
    PE_swing = r_swing.'*g_vec*m;
    PE_hip = r_hip.'*g_vec*m_h;
    PE=-PE_stance-PE_swing-PE_hip;
    
    %% M*ddq + C*dq + G = F + W^T*lambda
    % M
    
    M=simplify(jacobian(jacobian(KE,dq).',dq));
    
    % C
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
    
    % G
    G = simplify(jacobian(PE, q))';
    
    %% Constraints
    %
    h_t = x;
    h_n = y;
    h = [h_t; h_n];
    num_constraints = 2;
    
    W = jacobian(h, q);
    
    W_d = sym(zeros(num_constraints, N));
    for i = 1:num_constraints
        for j = 1:N
                W_d(i, j) = jacobian(W(i, j), q)*dq;
        end
    end
    
    A = simplify(W*inv(M)*W.');
    

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
fprintf(fid,'(q, q_d)\n');
fprintf(fid,'%% %s    Model of three-link robot cat.\n',...
	upper(fcn_name));

fprintf(fid,'%% Xianle Zeng\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[m, m_h, l, I_c, g, alpha, mu]=model_params();\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'x=q(1); y=q(2); theta_1=q(3); theta_2=q(4); \n');
fprintf(fid,'dx=q_d(1); dy=q_d(2); dtheta_1=q_d(3); dtheta_2=q_d(4); \n\n');

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
