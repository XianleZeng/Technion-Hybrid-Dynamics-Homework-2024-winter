clear

filename = strrep(mfilename,'gen_model','model_save');

if isempty(dir([filename,'.mat']))
    clear 

    filename = strrep(mfilename,'gen_model','model_save');

    %% variable declarations
    %
    syms th phi1 phi2 x y x1 y1 x2 y2 real
    syms th_d phi_d1 phi_d2 x_d y_d real
    syms th_dd phi_dd1 phi_dd2 x_dd y_dd real
    syms tau1 tau2
    syms m g l real
    
    %% generalized coordinates
    %
    qs = [phi1;phi2];
    qb = [x;y;th];  
    q = [qb; qs];
    b = length(qb);
    s = length(qs);
    
    %% first derivative of generalized coordinates
    %
    dqs = [th_d;phi_d1;phi_d2];
    dqb = [x_d; y_d];
    dq = [dqb; dqs];

    %% Second derivative of generalized coordinates
    %
    ddqs = [th_dd;phi_dd1;phi_dd2];
    ddqb = [x_dd; y_dd];
    ddq = [ddqb; ddqs];
    
    %% position of masses in system (center of mass for each link)
    %
    p_0 = [x; y];
    p_1 = p_0 - [l*(cos(phi1+th) + cos(th)); l*(sin(th) + sin(phi1+th))];
    p_2 = p_0 + [l*(cos(phi2+th) + cos(th)); l*(sin(th) + sin(phi2+th))];
    p_end_2 = p_0 + [l*cos(th); l*sin(th)];

    %% position of the center of mass 
    %
    p_cm = (1/(3*m))*(p_0*m + p_1*m + p_2*m);
    
    %% velocities of masses in system (center of mass for each link)
    %
    v_0 = jacobian(p_0,q)*dq;
    v_1 = jacobian(p_1,q)*dq;
    v_2 = jacobian(p_2,q)*dq;
    v_end_2 = jacobian(p_end_2,q)*dq;
    a_end_2 = jacobian(v_end_2,q)*ddq;

    %% velocity of the center of mass
    %
    v_cm = jacobian(p_cm,q)*dq;
    
    %% total angular momentum
    %
    I_c = m*l^2/3;
    temp_0 = cross([p_0 - p_cm; 0], [v_0 - v_cm; 0])*m;
    temp_1 = cross([p_1 - p_cm; 0], [v_1 - v_cm; 0])*m;
    temp_2 = cross([p_2 - p_cm; 0], [v_2 - v_cm; 0])*m;
    H_0 = temp_0(3) + I_c*th_d;
    H_1 = temp_1(3) + I_c*(th_d+phi_d1);
    H_2 = temp_2(3) + I_c*(th_d+phi_d2);
    H_total = H_0 + H_1 + H_2; 
    th_d_H = solve(H_total, th_d);

    %% kinetic energy of masses in system
    %
    KE_0 = simplify(m/2*(v_0'*v_0)) + (1/2)*I_c*th_d^2;
    KE_1 = simplify(m/2*(v_1'*v_1)) + (1/2)*I_c*(phi_d1+th_d)^2;
    KE_2 = simplify(m/2*(v_2'*v_2)) + (1/2)*I_c*(phi_d2+th_d)^2;
    
    %% total kinetic energy of system
    %
    KE = KE_1+KE_2+KE_0;
    
    %% potential energy of masses in system
    %
    PE_0 = p_0(2)*m*g;
    PE_1 = p_1(2)*m*g;
    PE_2 = p_2(2)*m*g;
    
    %% total potential energy of system
    %
    PE=PE_0+PE_1+PE_2;
    
    %%
    % M*ddq+C*dq+G = F
    M=simplify(jacobian(jacobian(KE,dq).',dq));
    
    %%
    N=max(size(q));
    syms C
    for k=1:N
	    for j=1:N
		    C(k,j)=0*g;
		    for i=1:N
			    C(k,j)=C(k,j)+1/2*(diff(M(k,j),q(i)) + diff(M(k,i),q(j)) - diff(M(i,j),q(k)))*dq(i);
		    end
	    end
    end
    B = C*dq;
    
    G=jacobian(PE,q).';
    
    %% Input Angle phi1, phi2
    syms t omega alfa Psi beta;
    
    t_f = 2*pi/omega;

    S = t^3*(6*t^2 - 15*t*t_f + 10*t_f^2)/t_f^5;
    phi1_input_1 = S*(-alfa+beta*sin(omega*t-Psi));
    phi2_input_1 = S*(alfa+beta*sin(omega*t+Psi));

    phi_input_1 = [phi1_input_1; phi2_input_1];
    phi_d_input_1 = diff(phi_input_1, t);
    phi_dd_input_1 = diff(phi_d_input_1, t);

    S = 1;
    phi1_input_2 = S*(-alfa+beta*sin(omega*t-Psi));
    phi2_input_2 = S*(alfa+beta*sin(omega*t+Psi));

    phi_input_2 = [phi1_input_2; phi2_input_2];
    phi_d_input_2 = diff(phi_input_2, t);
    phi_dd_input_2 = diff(phi_d_input_2, t);

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

N=max(size(q));

%% First, output model to a file called
%%
%%    dynamics_mat.m
%%
%

%% Output file header
%
fcn_name='dynamics_mat';
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function [M,B,G]=%s',fcn_name);
fprintf(fid,'(q, q_d)\n');
fprintf(fid,'%% %s    Model of three-link robot cat.\n',...
	upper(fcn_name));
fprintf(fid,'%% Input - values of generalized coordinate and their velocity q, q_dot \n');
fprintf(fid,'%% Output - matrices/vectors M(q), B(q, q_dot), G(q) for the dynamics equation of motion \n');

fprintf(fid,'%% Xianle Zeng\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[m, l, g]=model_params;\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'x=q(1); y=q(2); th=q(3); phi1=q(4); phi2=q(5);\n');
fprintf(fid,'x_d=q_d(1); y_d=q_d(2); th_d=q_d(3); phi_d1=q_d(4); phi_d2=q_d(5);\n\n');

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


%% Second, output model to a file called
%%
%%    angle_input.m
%%
%
%% Output file header
%
fcn_name='angles_input';
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function [phi, phi_d, phi_dd]=%s',fcn_name);
fprintf(fid,'(t)\n');
fprintf(fid,'%% %s    Calculate time profiles of input joint angles and their derivations\n',...
	upper(fcn_name));
fprintf(fid,'%% Input - time t \n');
fprintf(fid,'%% Output - vectors of phi(t), phi_d(t), phi_dd(t) \n');

fprintf(fid,'%% Xianle Zeng\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[omega,alfa,Psi,beta]=control_parameters;\n\n');
len_phi = max(size(phi_input_1));

%% Model output
%
fprintf(fid,'if (t <= 2*pi/omega) \n');
fprintf(fid,'\t%% phi \n');
fprintf(fid,'\tphi=zeros(%s,1);\n',num2str(2));
for k=1:len_phi
	if phi_input_1(k)~=0
		ttt=char(phi_input_1(k));
		fprintf(fid,'\tphi(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'\t%% phi_d \n');
fprintf(fid,'\tphi_d=zeros(%s,1);\n',num2str(2));
for k=1:len_phi
	if phi_d_input_1(k)~=0
		ttt=char(phi_d_input_1(k));
		fprintf(fid,'\tphi_d(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'\t%% phi_dd \n');
fprintf(fid,'\tphi_dd=zeros(%s,1);\n',num2str(2));
for k=1:len_phi
	if phi_dd_input_1(k)~=0
		ttt=char(phi_dd_input_1(k));
		fprintf(fid,'\tphi_dd(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'else \n');
fprintf(fid,'\t%% phi \n');
fprintf(fid,'\tphi=zeros(%s,1);\n',num2str(2));
for k=1:len_phi
	if phi_input_2(k)~=0
		ttt=char(phi_input_2(k));
		fprintf(fid,'\tphi(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'\t%% phi_d \n');
fprintf(fid,'\tphi_d=zeros(%s,1);\n',num2str(2));
for k=1:len_phi
	if phi_d_input_2(k)~=0
		ttt=char(phi_d_input_2(k));
		fprintf(fid,'\tphi_d(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'\t%% phi_dd \n');
fprintf(fid,'\tphi_dd=zeros(%s,1);\n',num2str(2));
for k=1:len_phi
	if phi_dd_input_2(k)~=0
		ttt=char(phi_dd_input_2(k));
		fprintf(fid,'\tphi_dd(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'end \n');

%% Third, output model to a file called
%%
%%    center_of_mass.m
%%
%
%% Output file header
%
fcn_name='center_of_mass';
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function [p_cm, v_cm]=%s',fcn_name);
fprintf(fid,'(q, q_d)\n');
fprintf(fid,'%% %s    center of mass of the cat\n',...
	upper(fcn_name));
fprintf(fid,'%% Input - values of generalized coordinate and their velocity q, q_dot \n');
fprintf(fid,'%% Output - position of the center of mass p_cm, velocity of the center of mass v_cm \n');

fprintf(fid,'%% Xianle Zeng\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[m, l, g]=model_params;\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'x=q(1); y=q(2); th=q(3); phi1=q(4); phi2=q(5);\n');
fprintf(fid,'x_d=q_d(1); y_d=q_d(2); th_d=q_d(3); phi_d1=q_d(4); phi_d2=q_d(5);\n\n');

%% Model output
%
fprintf(fid,'\n%% p_cm\n');
fprintf(fid, 'p_cm=zeros(2,1);\n');
for k=1:2
	if p_cm(k)~=0
		ttt=char(p_cm(k));
		fprintf(fid,'p_cm(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'\n%% v_cm\n');
fprintf(fid, 'v_cm=zeros(2,1);\n');
for k=1:2
	if v_cm(k)~=0
		ttt=char(v_cm(k));
		fprintf(fid,'v_cm(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'end \n');

%% Fourth, output model to a file called
%%
%%    Total_angular_momentum.m
%%
%
%% Output file header
%
fcn_name='Total_angular_momentum';
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function H_total=%s',fcn_name);
fprintf(fid,'(q, q_d, t)\n');
fprintf(fid,'%% %s    center of mass of the cat\n',...
	upper(fcn_name));
fprintf(fid,'%% Input - values of generalized coordinate and their velocity q, q_dot \n');
fprintf(fid,'%% Output - Total Angular Momentum of the System H_tatol \n');

fprintf(fid,'%% Xianle Zeng\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[m, l, g]=model_params;\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'x=q(1); y=q(2); th=q(3); phi1=q(4); phi2=q(5);\n');
fprintf(fid,'x_d=q_d(1); y_d=q_d(2); th_d=q_d(3); phi_d1=q_d(4); phi_d2=q_d(5);\n\n');
fprintf(fid,'qb = [x; y; th]; qb_d = [x_d; y_d; th_d];\n');

%% Get the tau
%
fprintf(fid, '[qb_dd, tau]=dyn_sol(qb,qb_d,t); \n');
fprintf(fid, 'tau1 = tau(1); tau2 = tau(2); \n');

%% Model output
%
fprintf(fid,'\n%% H_total\n');
fprintf(fid,'H_total=%s;\n',H_total);

fprintf(fid,'end \n');

%% Fifth, output model to a file called
%%
%%    theta_d_H.m
%%
%
%% Output file header
%
fcn_name='theta_d_H';
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function th_d_H=%s',fcn_name);
fprintf(fid,'(q, q_d)\n');
fprintf(fid,'%% %s    center of mass of the cat\n',...
	upper(fcn_name));

fprintf(fid,'%% Xianle Zeng\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[m, l, g]=model_params;\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'x=q(1); y=q(2); th=q(3); phi1=q(4); phi2=q(5);\n');
fprintf(fid,'x_d=q_d(1); y_d=q_d(2); th_d=q_d(3); phi_d1=q_d(4); phi_d2=q_d(5);\n\n');

%% Model output
%
fprintf(fid,'\n%% th_d_H\n');
fprintf(fid,'th_d_H=%s;\n',th_d_H);

fprintf(fid,'end \n');
