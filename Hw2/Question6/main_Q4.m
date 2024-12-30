clear;
close all;

%% Solve DAE 
%
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
tspan = 0:0.01:60;  
X0 = [0; 0; 0; 0; 0; 0];             
[t, X] = ode45(@state_eq_Q4, tspan, X0, options); 
[t_damping, X_damping] = ode45(@state_eq_Q4_damping, tspan, X0, options); 
n_data = length(t);
[m, l, I_c, d, w, b, ~]=model_params;

%% Input angles phi1 ,phi2 
%
phi_data = [];
phi_d_data = [];
phi_dd_data = [];
for i = 1:length(t)
    [phi, phi_d, phi_dd] = input_angle(t(i)); 
    phi_data = [phi_data, phi];
    phi_d_data = [phi_d_data, phi_d];
    phi_dd_data = [phi_dd_data, phi_dd];
end
phi = phi_data;
phi_d = phi_d_data;
phi_dd = phi_dd_data;

%% Reassign the variables (no damping)
%
x = X(:, 1)';
y = X(:, 2)';
theta = X(:, 3)';
x_d = X(:, 4)';
y_d = X(:, 5)';
theta_d = X(:, 6)';
q = [x; y; theta; phi];
q_d = [x_d; y_d; theta_d; phi_d];
qp = [x; y; theta];
qp_d = [x_d; y_d; theta_d];

%% Reassign the variables (damping)
%
x_damping = X_damping(:, 1)';
y_damping = X_damping(:, 2)';
theta_damping = X_damping(:, 3)';
x_damping_d = X_damping(:, 4)';
y_damping_d = X_damping(:, 5)';
theta_damping_d = X_damping(:, 6)';
q_damping = [x_damping; y_damping; theta_damping; phi];
q_damping_d = [x_damping_d; y_damping_d; theta_damping_d; phi_d];
qp_damping = [x_damping; y_damping; theta_damping];
qp_damping_d = [x_damping_d; y_damping_d; theta_damping_d];

%% Velocity of the point P projected along body direction e'_1 (no damping)
% 
P_d = [x_d(:), y_d(:)];
e1_prime = [cos(theta(:)), sin(theta(:))];
P_d_e1_prime = zeros(n_data, 1);
for i = 1:n_data
    P_d_e1_prime(i) = P_d(i,:)*e1_prime(i,:)';
end

%% Velocity of the point P projected along body direction e'_1 (damping)
% 
P_damping_d = [x_damping_d(:), y_damping_d(:)];
e1_prime_damping = [cos(theta_damping(:)), sin(theta_damping(:))];
P_d_e1_prime_damping = zeros(n_data, 1);
for i = 1:n_data
    P_d_e1_prime_damping(i) = P_damping_d(i,:)*e1_prime_damping(i,:)';
end

%% Lambda and qp_dd (damping)
%
lambda_damping_data = [];
qp_dd_data_damping = [];
damping = true;
for i = 1:length(t)
    [qp_dd_damping, F_qa_damping, lambda_damping]=dyn_sol_Q4(t(i),qp_damping(:,i)',qp_damping_d(:,i)',damping);
    lambda_damping_data = [lambda_damping_data, lambda_damping];
    qp_dd_data_damping = [qp_dd_data_damping, qp_dd_damping];
end
q_dd_damping = [qp_dd_data_damping; phi_dd];


%% Lambda and qp_dd (no damping)
%
lambda_data = [];
qp_dd_data = [];
damping = false;
for i = 1:length(t)
    [qp_dd, F_qa, lambda]=dyn_sol_Q4(t(i),qp(:,i)',qp_d(:,i)',damping);
    lambda_data = [lambda_data, lambda];
    qp_dd_data = [qp_dd_data, qp_dd];
end
q_dd = [qp_dd_data; phi_dd];

%% Center of mass (no damping)
%
r_cm_data = [];
v_cm_data = [];
a_cm_data = [];
for i = 1:length(t)
    [r_cm, v_cm, a_cm] = center_of_mass(q(:,i), q_d(:,i), q_dd(:, i));
    r_cm_data = [r_cm_data, r_cm];
    v_cm_data = [v_cm_data, v_cm];
    a_cm_data = [a_cm_data, a_cm];
end
x_cm = r_cm_data(1,:);
y_cm = r_cm_data(2,:);
vx_cm = v_cm_data(1,:);
vy_cm = v_cm_data(2,:);
ax_cm = a_cm_data(1,:);
ay_cm = a_cm_data(2,:);

%% Center of mass (damping)
%
r_cm_data_damping = [];
v_cm_data_damping = [];
a_cm_data_damping = [];
for i = 1:length(t)
    [r_cm_damping, v_cm_damping, a_cm_damping] = center_of_mass(q_damping(:,i), q_damping_d(:,i), q_dd_damping(:, i));
    r_cm_data_damping = [r_cm_data_damping, r_cm_damping];
    v_cm_data_damping = [v_cm_data_damping, v_cm_damping];
    a_cm_data_damping = [a_cm_data_damping, a_cm_damping];
end
x_cm_damping = r_cm_data_damping(1,:);
y_cm_damping = r_cm_data_damping(2,:);
vx_cm_damping = v_cm_data_damping(1,:);
vy_cm_damping = v_cm_data_damping(2,:);
ax_cm_damping = a_cm_data_damping(1,:);
ay_cm_damping = a_cm_data_damping(2,:);

%% inertia force (no damping)
%
inertia_force_x = zeros(n_data, 1);
for i = 1:n_data
    inertia_force_x(i) = ax_cm(1,i)*m;
end

%% inertia force (damping)
%
inertia_force_x_damping = zeros(n_data, 1);
for i = 1:n_data
    inertia_force_x_damping(i) = ax_cm_damping(1,i)*m;
end

%% reaction forces in x_hat dir  (no damp)
%
n_P = [cos(theta(:) + pi/2), sin(theta(:) + pi/2)];
n_3 = [cos(theta(:) + phi(:) + pi/2), sin(theta(:) + phi(:) + pi/2)];
t_3 = [cos(theta(:) + phi(:)), sin(theta(:) + phi(:))];
x_hat = [1, 0];

lambda_x = zeros(n_data,1);
for i = 1:n_data
    lambda_x(i) = lambda_data(1,i)*n_P(i,:)*x_hat' + lambda_data(2,i)*n_3(i,:)*x_hat'; 
end

reactant_x = lambda_x;


%% reaction forces in x_hat dir  (damping)
%
n_P_damping = [cos(theta_damping(:) + pi/2), sin(theta_damping(:) + pi/2)];
t_1_damping = [cos(theta_damping(:)), sin(theta_damping(:))];
t_2_damping = t_1_damping;
n_3_damping = [cos(theta_damping(:) + phi(:) + pi/2), sin(theta_damping(:) + phi(:) + pi/2)];
t_3_damping = [cos(theta_damping(:) + phi(:)), sin(theta_damping(:) + phi(:))];
x_hat = [1, 0];

lambda_x_damping = zeros(n_data,1);
for i = 1:n_data
    lambda_x_damping(i) = lambda_damping_data(1,i)*n_P_damping(i,:)*x_hat' + lambda_damping_data(2,i)*n_3_damping(i,:)*x_hat'; 
end

v_wheel_1_data_damping = zeros(2,n_data);
v_wheel_1_t_1_damping = zeros(n_data,1);

v_wheel_2_data_damping = zeros(2,n_data);
v_wheel_2_t_2_damping = zeros(n_data,1);

v_wheel_3_data_damping = zeros(2,n_data);
v_wheel_3_t_3_damping = zeros(n_data,1);

for i = 1:n_data
    [v_wheel_1_data_damping(:,i), v_wheel_2_data_damping(:,i), v_wheel_3_data_damping(:,i)] = v_wheel(q_damping(:,i), q_damping_d(:,i));
    v_wheel_1_t_1_damping(i) = t_1_damping(i,:)*v_wheel_1_data_damping(:,i);
    v_wheel_2_t_2_damping(i) = t_2_damping(i,:)*v_wheel_2_data_damping(:,i);
    v_wheel_3_t_3_damping(i) = t_3_damping(i,:)*v_wheel_3_data_damping(:,i);
end

F_damp_x_damping = zeros(n_data,1);

for i = 1:n_data
    F_damp_1_damping = -v_wheel_1_t_1_damping(i)*t_1_damping(i,:)';
    F_damp_2_damping = -v_wheel_2_t_2_damping(i)*t_2_damping(i,:)';
    F_damp_3_damping = -v_wheel_3_t_3_damping(i)*t_3_damping(i,:)';
    F_damp_x_damping(i) = F_damp_1_damping(1) +F_damp_2_damping(1) +  F_damp_3_damping(1);
end

reactant_x_damping = lambda_x_damping + F_damp_x_damping;


%% Plot Results %%%%%%%%%%%%%%%%%%%%%%%


%% Plot the angle phi
%
figure;
hold on; 
plot(t, rad2deg(phi), 'LineWidth', 4);
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('$\phi(t)$ [deg]', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 14);
% axis tight;
saveas(gcf, '../images/Q6_a.png');

%% Plot the angle theta
%
figure;
hold on;
plot(t, rad2deg(theta), 'LineWidth', 4, 'DisplayName', '$c = 0$');
hold on;
plot(t, rad2deg(theta_damping),'--', 'LineWidth', 4,  'DisplayName', '$c = 1$');
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('$\theta(t)$ [deg]', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 15;  
grid on;
grid on;
set(gca, 'FontSize', 14);
% axis tight;
saveas(gcf, '../images/Q6_b.png');

%% Plot velocity of the point P projected along body direction e'_1
%
figure;
hold on;
plot(t, P_d_e1_prime, 'LineWidth', 4, 'DisplayName', '$c = 0$');
hold on;
plot(t, P_d_e1_prime_damping, 'LineWidth', 4, 'DisplayName', '$c = 1$');
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('$v_P$ [m/s]', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 15;  
leg.Location='bestoutside';
grid on;
grid on;
set(gca, 'FontSize', 14);
% axis tight;
saveas(gcf, '../images/Q6_c.png');

%% Trajectory 
%
figure;
hold on;
plot(x, y, 'LineWidth', 4, 'DisplayName', '$c = 0$');
hold on;
plot(x_damping, y_damping, 'LineWidth', 4, 'DisplayName', '$c = 1$');
xlabel('x [m]', 'Interpreter', 'latex');
ylabel('y [m]', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 15;  
lgd.Location = 'southeast'; 
grid on;
grid on;
set(gca, 'FontSize', 14);
% axis tight;
saveas(gcf, '../images/Q6_d.png');


%% lambda
%
figure;
hold on;
plot(t, lambda_data(1,:), 'LineWidth', 3, 'DisplayName', '$\lambda_1$, c=0')
hold on;
plot(t, lambda_data(2,:), '--', 'LineWidth', 3, 'DisplayName', '$\lambda_2$, c=0')
hold on;
plot(t, lambda_damping_data(1,:), 'LineWidth', 3, 'DisplayName', '$\lambda_1$, c=1')
hold on;
plot(t, lambda_damping_data(2,:), '--', 'LineWidth', 3, 'DisplayName', '$\lambda_2$, c=1')
hold on
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('$\lambda$ [N]', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 15;  
grid on;
set(gca, 'FontSize', 14);
% axis tight;
saveas(gcf, '../images/Q6_e.png');


%% Inertia force and damping force
%
figure;
hold on;
plot(t, inertia_force_x, 'k', 'LineWidth', 3, 'DisplayName', 'Inertia force, c=0')
hold on;
plot(t, reactant_x, '--', 'Color', 'r', 'LineWidth', 3, 'DisplayName', 'Reactant, c=0')
hold on;
plot(t_damping, inertia_force_x_damping, 'Color', 'g', 'LineWidth', 3, 'DisplayName', 'Inertia force, c=1')
hold on;
plot(t_damping, reactant_x_damping,  '--', 'Color', 'b', 'LineWidth', 3, 'DisplayName', 'Reactant, c=1')
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('Force [N]', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 14;  
grid on;
set(gca, 'FontSize', 10);
saveas(gcf, '../images/Q6_f.png');

