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

    v_P = [x_d; y_d; 0];
    
    %% position (p_cm), velocities (v_cm) of center of mass
    %
    omega_theta = [0; 0; theta_d];
    omega_phi = [0; 0; phi_d];
    
    r_P_to_cm = [d*cos(theta); d*sin(theta); 0];
    r_cm = [x + d*cos(theta); y + d*sin(theta); 0];
    v_cm = [x_d; y_d; 0] + cross(omega_theta, r_P_to_cm);
    
    %% Rotation matrix of the float-based coordinate:
    %
    R = [cos(theta), -sin(theta), 0;
           sin(theta), cos(theta), 0;
           0, 0, 1];
    
    %% position, velocities of the wheels
    %
    r_cm_to_wheel_1 = [-d; w/2; 0];
    r_cm_to_wheel_2 = [-d; -w/2; 0];
    r_cm_to_wheel_3 = [l-d+b*cos(phi); b*sin(phi); 0];
    
    r_wheel_1 =  R*r_cm_to_wheel_1 + r_cm;
    v_wheel_1 = v_cm + cross(omega_theta, r_cm_to_wheel_1);
    
    r_wheel_2 =  R*r_cm_to_wheel_2 + r_cm;
    v_wheel_2 = v_cm + cross(omega_theta, r_cm_to_wheel_2);
    
    r_wheel_3 =  R*r_cm_to_wheel_3 + r_cm;
    v_wheel_3 = v_cm + cross(omega_theta + omega_phi, r_cm_to_wheel_3);
    %% kinetic energy and ﻿dissipation function
    %
    KE = simplify(m/2*(v_cm'*v_cm) + (1/2)*I_c*theta_d^2);
    D = (1/2)*c*norm(v_wheel_1)^2 + (1/2)*c*norm(v_wheel_2)^2 + (1/2)*c*norm(v_wheel_3)^2;
    
    %% Nonholonomic Constraints 
    %
    e2_prime = [0; 1; 0];
    n1_hat = R*e2_prime;
    n3_hat = - [1; 0; 0];
    
    f_1 = simplify(v_P'*n1_hat);
    f_3 = simplify(v_wheel_3'*n3_hat);
    
    f_1 = collect(collect(collect(collect(f_1, x_d), y_d),theta_d), phi_d);
    f_3 = collect(collect(collect(collect(f_3, x_d), y_d),theta_d), phi_d);
    
    f = [f_1; f_3];
    m = max(size(f));
   

    W = sym(zeros(m, N));
    for i = 1:m
        for j = 1:N
            [coeff, ~] = coeffs(f(i), dq(j));
            if (length(coeff) == 2)
                W(i, j) = coeff(1);
            end
        end
    end

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
fprintf(fid,'function [M,Mpp,Mpa,Maa,B,Bp,Ba,W,W_d,Wp,Wa,Wp_d,Wa_d]=%s',fcn_name);
fprintf(fid,'(q, q_d)\n');
fprintf(fid,'%% %s    Model of three-link robot cat.\n',...
	upper(fcn_name));

fprintf(fid,'%% Xianle Zeng\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[m, l, I_c, d, w, b, c]=model_params;\n\n');

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