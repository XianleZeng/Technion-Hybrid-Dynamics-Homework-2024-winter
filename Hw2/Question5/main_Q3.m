clear;
close all;

%% Solve DAE without damping
%
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
tspan = 0:0.01:100;  
X0 = [0; 0; 0; 0; 0; 0; 0; 0];   % x_0, y_0, theta_0, phi_0, dx_0, dy_0, dtheta_0, dphi_0     
[t, X] = ode45(@state_eq_Q3, tspan, X0, options);  % dX = f(X)
n_data = length(t);
damping = false;
[m, l, I_c, d, w, b, ~]=model_params;

%% Reassign the variables
%
x = X(:, 1)';
y = X(:, 2)';
theta = X(:, 3)';
phi = X(:, 4)';
x_d = X(:, 5)';
y_d = X(:, 6)';
theta_d = X(:, 7)';
phi_d = X(:, 8)';
q = [x; y; theta; phi];
q_d = [x_d; y_d; theta_d; phi_d];
qp = [x; y; theta];
qp_d = [x_d; y_d; theta_d];

%% Velocity of the point P projected along body direction e'_1
%
P_d = [x_d(:), y_d(:)];
e1_prime = [cos(theta(:)), sin(theta(:))];
P_d_e1_prime = zeros(n_data, 1);
for i = 1:n_data
    P_d_e1_prime(i) = P_d(i,:)*e1_prime(i,:)';
end

%% Lambda and qp_dd
%
lambda_data = [];
q_dd_data = [];
for i = 1:length(t)
    [q_dd, lambda]=dyn_sol_Q3(t(i),q(:,i),q_d(:,i));
    lambda_data = [lambda_data, lambda];
    q_dd_data = [q_dd_data, q_dd];
end
q_dd = q_dd_data;

%% Center of mass
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

%% inertia force
%
inertia_force_x = zeros(n_data, 1);
for i = 1:n_data
    inertia_force_x(i) = ax_cm(1,i)*m;
end

%% reaction forces in x_hat dir
%
n_P = [cos(theta(:) + pi/2), sin(theta(:) + pi/2)];
n_3 = [cos(theta(:) + phi(:) + pi/2), sin(theta(:) + phi(:) + pi/2)];
x_hat = [1, 0];

lambda_x = zeros(n_data,1);
for i = 1:n_data
    lambda_x(i) = lambda_data(1,i)*n_P(i,:)*x_hat' + lambda_data(2,i)*n_3(i,:)*x_hat'; 
end

F_damp_x = 0;

reactant_x = lambda_x + F_damp_x;


%%%%%%%%%%%%% plot
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
saveas(gcf, 'images/Q3_a.png');

%% Plot the angle theta
%
figure;
hold on;
plot(t, rad2deg(theta), 'LineWidth', 4);
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('$\theta(t)$ [deg]', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 14);
% axis tight;
saveas(gcf, 'images/Q3_b.png');

%% Plot velocity of the point P projected along body direction e'_1
%
figure;
hold on;
plot(t, P_d_e1_prime, 'LineWidth', 4);
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('$v_P$ [m/s]', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 14);
% axis tight;
saveas(gcf, 'images/Q3_c.png');

%% trajectory
%
figure;
hold on;
plot(x, y, 'LineWidth', 4);
xlabel('x [m]', 'Interpreter', 'latex');
ylabel('y [m]', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 14);
% axis tight;
saveas(gcf, 'images/Q3_d.png');


%% lambda
%
figure;
hold on;
plot(t, lambda_data(1,:), 'LineWidth', 3, 'DisplayName', '$\lambda_1$')
hold on;
plot(t, lambda_data(2,:), 'LineWidth', 3, 'DisplayName', '$\lambda_2$')
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('$\lambda$ [N]', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 15;  
grid on;
set(gca, 'FontSize', 14);
% axis tight;
saveas(gcf, 'images/Q3_e.png');


%% Inertia force and damping force
%
figure;
hold on;
plot(t, inertia_force_x, 'LineWidth', 3, 'DisplayName', 'Inertia force')
hold on;
plot(t, reactant_x, '--', 'LineWidth', 3, 'DisplayName', 'Reactant')
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('Forces [N]', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 15;  
grid on;
set(gca, 'FontSize', 14);
saveas(gcf, 'images/Q4_f.png');


%% slip velocity
%
v_slip_1 = zeros(length(t), 1);
v_slip_2 = zeros(length(t), 1);

for i = 1:length(t)
    damping = false;
    [M,Mpp,Mpa,Maa,B,Bp,Ba,G,Gp,Ga,W,W_d,Wp,Wa,Wp_d,Wa_d]=dynamics_mat(q(:, i), q_d(:, i), damping);
    w1 = W(1, :);
    w2 = W(2, :);
    v_slip_1(i) = w1 * q_d(:, i);
    v_slip_2(i) = w2 * q_d(:, i);
end

figure;
plot(t, v_slip_1, 'LineWidth', 3);hold on;
plot(t, v_slip_2, '--', 'LineWidth', 2); 
ylim([-2, 2]);
grid on;
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Slip Velocities $[\mathrm{m/s}]$', 'Interpreter', 'latex', 'FontSize', 15);
legend({'$\mathbf{w}_1(\mathbf{q}(t)) \dot{\mathbf{q}}(t)$', ...
    '$\mathbf{w}_2(\mathbf{q}(t)) \dot{\mathbf{q}}(t)$'}, ...
    'Interpreter', 'latex', 'FontSize', 15, 'Location', 'northeast');
set(gca, 'FontSize', 15, 'LineWidth', 1.2);
ylim([-2, 2]);
saveas(gcf, 'images/Q4_g.png');