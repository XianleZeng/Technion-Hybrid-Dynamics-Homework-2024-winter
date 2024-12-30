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
    syms m l I_c real
    syms d w b real
    syms c real
    
    %% generalized coordinates
    %
    q_p = [x; y; theta];
    q_a = phi;
    dq_p = [x_d; y_d; theta_d];
    dq_a = phi_d;
    ddq_p = [x_dd; y_dd; theta_dd];
    ddq_a = phi_dd;

    q = [q_p; q_a];
    dq = [dq_p; dq_a];
    ddq = [ddq_p; ddq_a];

    N=max(size(q));
    n_p = length(q_p);
    n_a = length(q_a);

    r_P = [x; y; 0];
    v_P = [x_d; y_d; 0];
    
    %% position (p_cm), velocities (v_cm) of center of mass
    %
    omega_theta = [0; 0; theta_d];
    omega_phi = [0; 0; phi_d];
    
    r_cm = [x + d*cos(theta); y + d*sin(theta)];
    v_cm = jacobian(r_cm, q)*dq;
    a_cm = jacobian(v_cm, [q; dq])*[dq; ddq]; 
    
    %% Rotation matrix of the float-based coordinate:
    %
    R = [cos(theta), -sin(theta), 0;
           sin(theta), cos(theta), 0;
           0, 0, 1];
    R_phi = [cos(phi+theta), -sin(phi+theta), 0;
           sin(phi+theta), cos(phi+theta), 0;
           0, 0, 1];
    
    %% position, velocities of the wheels
    %
    r_P_to_wheel_1 = [w/2*cos(theta+pi/2); w/2*sin(theta+pi/2); 0];
    r_P_to_wheel_2 = [-w/2*cos(theta+pi/2); -w/2*sin(theta+pi/2); 0];
    r_P_to_wheel_3 = [l*cos(theta)+b*cos(theta+phi); l*sin(theta)+b*sin(theta+phi); 0];
    
    r_wheel_1 =  r_P_to_wheel_1 + r_P;
    v_wheel_1 = jacobian(r_wheel_1, q)*dq;
    
    r_wheel_2 =  r_P_to_wheel_2 + r_P;
    v_wheel_2 = jacobian(r_wheel_2, q)*dq;
    
    r_wheel_3 =  r_P_to_wheel_3 + r_P;
    v_wheel_3 = jacobian(r_wheel_3, q)*dq;


    %% kinetic energy and ﻿dissipation function
    %
    KE = simplify(m/2*(v_cm'*v_cm) + (1/2)*I_c*theta_d^2);
    D = (1/2)*c*norm(v_wheel_1)^2 + (1/2)*c*norm(v_wheel_2)^2 + (1/2)*c*norm(v_wheel_3)^2;
    
    %% Nonholonomic Constraints 
    %
    e2_prime = [0; 1; 0];
    n1_hat = R*e2_prime;
    n3_hat = R_phi*[0; 1; 0];
    
    f_1 = simplify(v_P'*n1_hat);
    f_3 = simplify(v_wheel_3'*n3_hat);
    
    f = [f_1; f_3];
    m = max(size(f));
    
    W = jacobian(f, dq);

    Wp = W(:, 1:n_p);
    Wa = W(:, n_p+1:end);

    W_d = sym(zeros(m, N));
    for i = 1:m
        for j = 1:N
                W_d(i, j) = jacobian(W(i, j), q)*dq;
        end
    end

    Wp_d = W_d(:, 1:n_p);
    Wa_d = W_d(:, n_p+1:end);
    
    %%
    % M*ddq+B = F
    M=simplify(jacobian(jacobian(KE,dq).',dq));
    Mpp = M(1:n_p, 1:n_p);
    Mpa = M(1:n_p, n_p+1:end);
    Maa = M(n_p+1:end, n_p+1:end);

    %%
    syms C1
    for k=1:N
        for j=1:N
	        C1(k,j)=0;
	        for i=1:N
		        C1(k,j)=C1(k,j)+1/2*(diff(M(k,j),q(i)) + diff(M(k,i),q(j)) - diff(M(i,j),q(k)))*dq(i);
	        end
        end
    end
    C2 = simplify(jacobian(D, dq).');
    B = C1*dq + C2;

    Bp = B(1:n_p);
    Ba = B(n_p+1:end);

    syms G
    for j=1:N
        G(j) = 0;
    end
    G = G';

    Ga = G(1:n_a);
    Gp = G(n_a+1:end);

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
fprintf(fid,'function [M,Mpp,Mpa,Maa,B,Bp,Ba,G,Gp,Ga,W,W_d,Wp,Wa,Wp_d,Wa_d]=%s',fcn_name);
fprintf(fid,'(q, q_d, damping)\n');
fprintf(fid,'%% %s    Model of three-link robot cat.\n',...
	upper(fcn_name));

fprintf(fid,'%% Xianle Zeng\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[m, l, I_c, d, w, b, c]=model_params(''damping'', damping);\n\n');

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

fprintf(fid,'%% Mpp matrix\n');
fprintf(fid,'Mpp=zeros(%s);\n',num2str(n_p));
for k=1:n_p
	for j=1:n_p
		if Mpp(k,j)~=0
			ttt=char(Mpp(k,j));
			fprintf(fid,'Mpp(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'%% Maa matrix\n');
fprintf(fid,'Maa=zeros(%s);\n',num2str(n_a));
for k=1:n_a
	for j=1:n_a
		if Maa(k,j)~=0
			ttt=char(Maa(k,j));
			fprintf(fid,'Maa(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'%% Mpa matrix\n');
fprintf(fid,'Mpa=zeros(%s, %s);\n',num2str(n_p), num2str(n_a));
for k=1:n_p
	for j=1:n_a
		if Mpa(k,j)~=0
			ttt=char(Mpa(k,j));
			fprintf(fid,'Mpa(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
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

fprintf(fid,'\n%% Bp matrix\n');
fprintf(fid,'Bp=zeros(%s,1);\n',num2str(n_p));
for k=1:n_p
	if Bp(k)~=0
		ttt=char(Bp(k));
		fprintf(fid,'Bp(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'\n%% Ba matrix\n');
fprintf(fid,'Ba=zeros(%s,1);\n',num2str(n_a));
for k=1:n_a
	if Ba(k)~=0
		ttt=char(Ba(k));
		fprintf(fid,'Ba(%s)=%s;\n',num2str(k),ttt);
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

fprintf(fid,'\n%% Ga matrix\n');
fprintf(fid,'Ga=zeros(%s,1);\n',num2str(n_a));
for k=1:n_a
	if Ga(k)~=0
		ttt=char(Ga(k));
		fprintf(fid,'Ga(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'\n%% Gp matrix\n');
fprintf(fid,'Gp=zeros(%s,1);\n',num2str(n_p));
for k=1:n_p
	if Gp(k)~=0
		ttt=char(Gp(k));
		fprintf(fid,'Gp(%s)=%s;\n',num2str(k),ttt);
	end
end


fprintf(fid,'\n%% W matrix\n');
fprintf(fid,'W=zeros(%s, %s);\n',num2str(m), num2str(N));
for k = 1:m
    for j = 1:N
		if W(k,j)~=0
			ttt=char(W(k,j));
			fprintf(fid,'W(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
    end
end

fprintf(fid,'\n%% Wp matrix\n');
fprintf(fid,'Wp=zeros(%s, %s);\n',num2str(m), num2str(n_p));
for k = 1:m
    for j = 1:n_p
		if Wp(k,j)~=0
			ttt=char(Wp(k,j));
			fprintf(fid,'Wp(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
    end
end

fprintf(fid,'\n%% Wa matrix\n');
fprintf(fid,'Wa=zeros(%s,1);\n',num2str(m));
for k=1:m
	if Wa(k)~=0
		ttt=char(Wa(k));
		fprintf(fid,'Wa(%s)=%s;\n',num2str(k),ttt);
	end
end


fprintf(fid,'\n%% W_d matrix\n');
fprintf(fid,'W_d=zeros(%s, %s);\n',num2str(m), num2str(N));
for k = 1:m
    for j = 1:N
		if W_d(k,j)~=0
			ttt=char(W_d(k,j));
			fprintf(fid,'W_d(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
    end
end

fprintf(fid,'\n%% Wp_d matrix\n');
fprintf(fid,'Wp_d=zeros(%s, %s);\n',num2str(m), num2str(n_p));
for k = 1:m
    for j = 1:n_p
		if Wp_d(k,j)~=0
			ttt=char(Wp_d(k,j));
			fprintf(fid,'Wp_d(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
    end
end

fprintf(fid,'\n%% Wa_d matrix\n');
fprintf(fid,'Wa_d=zeros(%s,1);\n',num2str(m));
for k=1:m
	if Wa_d(k)~=0
		ttt=char(Wa_d(k));
		fprintf(fid,'Wa_d(%s)=%s;\n',num2str(k),ttt);
	end
end

%% Second, output model to a file called
%%
%%    center_of_mass.m
%%
%
%% Output file header
%
fcn_name='center_of_mass';
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function [r_cm, v_cm, a_cm]=%s',fcn_name);
fprintf(fid,'(q, q_d, q_dd)\n');
fprintf(fid,'%% %s    center of mass of the cat\n',...
	upper(fcn_name));
fprintf(fid,'%% Input - values of generalized coordinate and their velocity q, q_dot \n');
fprintf(fid,'%% Output - position of the center of mass p_cm, velocity of the center of mass v_cm \n');

fprintf(fid,'%% Xianle Zeng\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[m, l, I_c, d, w, b, c]=model_params;\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'x=q(1); y=q(2); theta=q(3); phi=q(4);\n');
fprintf(fid,'x_d=q_d(1); y_d=q_d(2); theta_d=q_d(3); phi_d=q_d(4);\n\n');
fprintf(fid,'x_dd=q_dd(1); y_dd=q_dd(2); theta_dd=q_dd(3); phi_dd=q_dd(4);\n\n');

%% Model output
%
fprintf(fid,'\n%% r_cm\n');
fprintf(fid, 'r_cm=zeros(2,1);\n');
for k=1:2
	if r_cm(k)~=0
		ttt=char(r_cm(k));
		fprintf(fid,'r_cm(%s)=%s;\n',num2str(k),ttt);
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

fprintf(fid,'\n%% a_cm\n');
fprintf(fid, 'a_cm=zeros(2,1);\n');
for k=1:2
	if a_cm(k)~=0
		ttt=char(a_cm(k));
		fprintf(fid,'a_cm(%s)=%s;\n',num2str(k),ttt);
	end
end


fprintf(fid,'end \n');




%% Third, output model to a file called
%%
%%    v_wheel_3.m
%%
%
%% Output file header
%
fcn_name='v_wheel_3';
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function v_wheel_3=%s',fcn_name);
fprintf(fid,'(q, q_d)\n');
fprintf(fid,'%% Xianle Zeng\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[m, l, I_c, d, w, b, c]=model_params;\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'x=q(1); y=q(2); theta=q(3); phi=q(4);\n');
fprintf(fid,'x_d=q_d(1); y_d=q_d(2); theta_d=q_d(3); phi_d=q_d(4);\n\n');

%% Model output
%
fprintf(fid,'\n%% v_wheel_3\n');
fprintf(fid, 'v_wheel_3=zeros(2,1);\n');
for k=1:2
	if v_wheel_3(k)~=0
		ttt=char(v_wheel_3(k));
		fprintf(fid,'v_wheel_3(%s)=%s;\n',num2str(k),ttt);
	end
end
