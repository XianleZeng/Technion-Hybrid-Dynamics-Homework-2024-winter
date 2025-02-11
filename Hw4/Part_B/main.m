close all
clc 
clear
 
save = 1;

[m_1, m_2, l, h, J_1, J_2, R, g, mu]=model_params();
phi_d_0 = 17.45*0.25 + 1.291*0.75; % rad/s
q0 = [0; R; 0; 0; 0; 0; 0; phi_d_0];
t0 = 0;
 
time_span = 10; % sec
 
t_total = [];
q_total = [];
 
contact = 1; % assume begin with contact
slip = 0;
separation = 0;

contact_type = [];
t_contact_type = [0];

Ie_stick = -1;
Ie_slip = -1;
 
while (separation == 0 & ~isempty(Ie_stick) & ~isempty(Ie_slip))
    if (contact == 1)
        op_stick = odeset('reltol',1e-8,'abstol',1e-8,'Events',@(t, q)events_stick(t, q));
        [t_stick, q_stick, Te_stick, Xe_stick, Ie_stick] = ode45(@(t,q)sys_stick(t, q), [t0, time_span], q0, op_stick);
        
        t_total = [t_total; t_stick];
        q_total = [q_total; q_stick];
        
        t0 = t_stick(end);
        q0 = Xe_stick;
        
        t_contact_type = [t_contact_type; Te_stick];
        contact_type = [contact_type; 1]; % STICK
        
        if ~isempty(Ie_stick)
           contact = 0;
           slip = 1;
           slip_dir = Ie_stick;
        end
    end


    if (slip == 1)
        op_slip = odeset('reltol',1e-8,'abstol',1e-8,'Events',@(t, q)events_slip(t, q, slip_dir));
        [t_slip, q_slip, Te_slip, Xe_slip, Ie_slip] = ode45(@(t, q)sys_slip(t, q, slip_dir), [t0, time_span], q0, op_slip);
        if (length(Ie_slip) == 2)
            Ie_slip = Ie_slip(2);
            Te_slip = Te_slip(2);
            Xe_slip = Xe_slip(2, :);
        end

        t_total = [t_total; t_slip];
        q_total = [q_total; q_slip];
        
        t0 = t_slip(end);
        q0 = Xe_slip';
        
        t_contact_type = [t_contact_type; Te_slip];
    
    
        if (slip_dir == 1)
            contact_type = [contact_type; 2]; % SLIP (positive)
        elseif (slip_dir == 2)
            contact_type = [contact_type; 3]; % SLIP (negative)
        end
    
        if (Ie_slip == 1)
        
            q_ = q_slip(end, 1:4)';  q_d_ = q_slip(end, 5:8)';
            
            [~, lambda] = dyn_sol_stick(t_slip(end),q_,q_d_);  % A*[q_dd; lambda] = B
            
            lambda_t = lambda(1);
            lambda_n = lambda(2);
        
            if (abs(lambda_t) < mu*lambda_n)
                contact = 1;
                slip = 0;
            else
                if slip_dir == 1
                    slip_dir = 2;
                elseif slip_dir == 2
                    slip_dir = 1;
                end  
            end
        
        elseif (Ie_slip == 2)
            separation = 1
            t_contact_type = [t_contact_type; Te_slip];
            contact_type = [contact_type; 4]; % STICK
        end
    end 
end

%%
t = t_total;

q = q_total(:, 1:4)';  q_d = q_total(:, 5:8)';
x = q(1,:);
y = q(2,:);
th = q(3,:);
ph = q(4,:);
x_d = q_d(1,:);
y_d = q_d(2,:);
th_d = q_d(3,:);
ph_d = q_d(4,:);


%% plot (a): theta(t) and phi(t)
figure(1)
 
subplot(2,1,1)
% plot(t, rad2deg(th), '-','LineWidth', 2)

t_con = [t_contact_type; time_span];
for i = 1:length(contact_type)
    
    con = contact_type(i);
    
    idx = find(t >= t_con(i)-0.0001 & t < t_con(i+1));
    if isempty(idx)
        break;
    end
    if (idx(1)~=1)
        idx = [idx(1)-1; idx];
    end
    t_i = t(idx);
    th_i = rad2deg(th(idx)); 
    
    if con == 1
        plot(t_i, th_i, '-', 'Color', '#0072BD', 'LineWidth', 1.5)
    elseif con == 2 || con == 3
        plot(t_i, th_i, '--', 'Color', '#EDB120', 'LineWidth', 1.5)
    end
    
    hold on
end

legend("${\theta(t)}$ - Stick", "${\theta(t)}$ - Slip",'Interpreter','Latex')
 
xlabel('${t}$ [sec]','Interpreter','Latex')
ylabel('${\theta}$ [deg]','Interpreter','Latex')
title('${\theta(t)}$','Interpreter','Latex')
 
subplot(2,1,2)
% plot(t, rad2deg(ph), '-','LineWidth', 2)

t_con = [t_contact_type; time_span];
for i = 1:length(contact_type)
    
    con = contact_type(i);
    
    idx = find(t >= t_con(i)-0.0001 & t < t_con(i+1));
    if isempty(idx)
        break;
    end
    if (idx(1)~=1)
        idx = [idx(1)-1; idx];
    end
    t_i = t(idx);
    ph_i = rad2deg(ph(idx)); 
    
    if con == 1
        plot(t_i, ph_i, '-', 'Color', '#0072BD', 'LineWidth', 1.5)
    elseif con == 2 || con == 3
        plot(t_i, ph_i, '--', 'Color', '#EDB120', 'LineWidth', 1.5)
    end
    
    hold on
end

legend("${\phi(t)}$ - Stick", "${\phi(t)}$ - Slip",'Interpreter','Latex')

xlabel('${t}$ [sec]','Interpreter','Latex')
ylabel('${\phi}$ [deg]','Interpreter','Latex')
title('${\phi(t)}$','Interpreter','Latex')
if (save == 1)
    saveas(gcf, ".\images\a_omega_0_"+ num2str(phi_d_0) + ".png");
end

%% plot (b): vt(t)
figure(2)
 
vt = x_d-R*ph_d;
plot(t, vt, '-','LineWidth', 2)
 
xlabel('${t}$ [sec]','Interpreter','Latex')
ylabel('${v_t}$ [m/s]','Interpreter','Latex')
title('${v_t(t)}$','Interpreter','Latex')
if (save == 1)
    saveas(gcf, ".\images\b_omega_0_"+ num2str(phi_d_0) + ".png");
end
%% plot (c): lam_n(t)
figure(3)
[lam_n, ~] = lambda_cal(t, q, q_d, t_contact_type, contact_type, separation);
plot(t, lam_n, '-','LineWidth', 2)
t_con = [t_contact_type; time_span];
for i = 1:length(contact_type)
    
    con = contact_type(i);
    
    idx = find(t >= t_con(i)-0.0001 & t < t_con(i+1));
    if isempty(idx)
        break;
    end
    if (idx(1)~=1)
        idx = [idx(1)-1; idx];
    end
    t_i = t(idx);
    lam_n_i = lam_n(idx); 
    
    if con == 1
        plot(t_i, lam_n_i, '-', 'Color', '#0072BD', 'LineWidth', 1.5)
    elseif con == 2 || con == 3
        plot(t_i, lam_n_i, '--', 'Color', '#EDB120', 'LineWidth', 1.5)
    end
    
    hold on
end

legend("${\lambda_n(t)}$ - Stick", "${\lambda_n(t)}$ - Slip",'Interpreter','Latex')


hold on
plot(t, t*0, '--', 'Color', '#A2142F', 'LineWidth', 1.5)

xlabel('${t}$ [sec]','Interpreter','Latex')
ylabel('${\lambda_n}$ [N]','Interpreter','Latex')
title('${\lambda_n(t)}$','Interpreter','Latex')
if (save == 1)
    saveas(gcf, ".\images\c_omega_0_"+ num2str(phi_d_0) + ".png");
end

%% plot (d): lam_t(t)/lam_n(t)
figure(4)
[lam_n, lam_t] = lambda_cal(t, q, q_d, t_contact_type, contact_type, separation);
plot(t, lam_t./lam_n, '-', 'Color', '#0072BD', 'LineWidth', 2)
hold on
plot(t, mu*ones(length(t),1), '--', 'Color', '#A2142F', 'LineWidth', 1.5)
plot(t, -mu*ones(length(t),1), '--', 'Color', '#A2142F', 'LineWidth', 1.5)
 
xlabel('${t}$ [sec]','Interpreter','Latex')
ylabel('${\lambda_t/\lambda_n}$ [-]','Interpreter','Latex')
title('${\lambda_t(t)/\lambda_n(t)}$','Interpreter','Latex')
 
legend("${\lambda_t(t)/\lambda_n(t)}$", "${\pm \mu}$",'Interpreter','Latex')
if (save == 1)
    saveas(gcf, ".\images\d_omega_0_"+ num2str(phi_d_0) + ".png");
end

%% plot (e): x(t)
figure(5)
 
t_con = [t_contact_type; time_span];
for i = 1:length(contact_type)
    
    con = contact_type(i);
    
    idx = find(t >= t_con(i)-0.0001 & t < t_con(i+1));
    if isempty(idx)
        break;
    end
    if (idx(1)~=1)
        idx = [idx(1)-1; idx];
    end
    t_i = t(idx);
    x_i = x(idx); 
    
    if con == 1
        plot(t_i, x_i, '-', 'Color', '#0072BD', 'LineWidth', 1.5)
    elseif con == 2 || con == 3
        plot(t_i, x_i, '--', 'Color', '#EDB120', 'LineWidth', 1.5)
    end
    
    hold on
end
 
xlabel('${t}$ [sec]','Interpreter','Latex')
ylabel('${x}$ [m]','Interpreter','Latex')
title('${x(t)}$','Interpreter','Latex')
 
legend("${x(t)}$ - Stick", "${x(t)}$ - Slip",'Interpreter','Latex')
if (save == 1)
    saveas(gcf, ".\images\e_omega_0_"+ num2str(phi_d_0) + ".png");
end