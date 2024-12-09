clear;
close all;

%% Solve DAE
%
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
tspan = 0:0.01:40;  
X0 = [0; 0; 0; 0; 0; 0];             
[t, X] = ode45(@state_eq, tspan, X0, options); 

%% Input angles phi1 ,phi2
%
phi_data = [];
phi_d_data = [];
for i = 1:length(t)
    [phi, phi_d, ~] = angles_input(t(i)); 
    phi_data = [phi_data, phi];
    phi_d_data = [phi_d_data, phi_d];
end
phi1 = phi_data(1,:)';
phi2 = phi_data(2,:)';
phi1_d = phi_d_data(1,:)';
phi2_d = phi_d_data(2,:)';

%% Reassign the variables
%
x = X(:, 1);
y = X(:, 2);
theta = X(:, 3);
x_d = X(:, 4);
y_d = X(:, 5);
theta_d = X(:, 6);
q = [x, y, theta, phi1, phi2];
q_d = [x_d, y_d, theta_d, phi1_d, phi2_d];
qb = [x, y, theta];
qb_d = [x_d, y_d, theta_d];
w = 1; %% Not sure, need to check!

%% Center of mass
%
p_cm_data = [];
v_cm_data = [];
for i = 1:length(t)
    [p_cm, v_cm] = center_of_mass(q(i,:), q_d(i,:));
    p_cm_data = [p_cm_data, p_cm];
    v_cm_data = [v_cm_data, v_cm];
end
x_cm = p_cm_data(1,:);
y_cm = p_cm_data(2,:);


%% Torque and qb_dd
%
tau_data = [];
qb_dd_data = [];
for i = 1:length(t)
    [qb_dd, tau]=dyn_sol(qb(i,:)',qb_d(i,:)',t(i));
    tau_data = [tau_data, tau];
    qb_dd_data = [qb_dd_data, qb_dd];
end
theta_dd = qb_dd_data(3,:)';

%%


%% Plot the angle of middle link 
%
figure;
hold on;
plot(t, theta, 'LineWidth', 2);
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('$\theta(t)$ [rad]', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 14);
axis tight;
saveas(gcf, 'images/a.png');

%% Figure for horizontal position of the middle link's center x(t)/l 
% With center-of-mass position x_cm(t)/l
% Need to plot: H(t)/ml^2w
%
figure;
hold on;
l = 0.1;
plot(t, x/l, 'LineWidth', 2, 'DisplayName', '$\frac{x(t)}{l}$');
hold on
plot(t, x_cm/l, '--', 'LineWidth', 2, 'DisplayName', '$\frac{x_c(t)}{l}$');
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('Normalized horizontal position [m/m]', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
grid on;
set(gca, 'FontSize', 14);
axis tight;
saveas(gcf, 'images/b.png');

%% Figure for horizontal velocity of the middle link's center x_dot(t)/lw 
% With center-of-mass position x_dot_cm(t)/lw
%
figure;
hold on;
l = 0.1;
plot(t, x_d/(l*w), 'LineWidth', 2, 'DisplayName', '$\frac{\dot{x}(t)}{l \omega}$');
hold on
plot(t, v_cm_data(1,:)'/(l*w), '--', 'LineWidth', 2, 'DisplayName', '$\frac{\dot{x}_c(t)}{l \omega}$');
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('Normalized horizontal velocity $[\frac{m/s}{m/s}]$', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
grid on;
set(gca, 'FontSize', 14);
axis tight;
saveas(gcf, 'images/c.png');

%% Figure for vertical velocity of the middle link's center y_dot(t)/lw 
% With center-of-mass position y_dot_cm(t)/lw
%
figure;
hold on;
l = 0.1;
plot(t, y_d/l, 'LineWidth', 2, 'DisplayName', '$\frac{\dot{y}(t)}{l \omega}$');
hold on
plot(t, v_cm_data(2,:)'/(w*l), '--', 'LineWidth', 2, 'DisplayName', '$\frac{\dot{y}_c(t)}{l \omega}$');
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('Normalized vertical velocity $[\frac{m/s}{m/s}]$', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
grid on;
set(gca, 'FontSize', 14);
axis tight;
saveas(gcf, 'images/d.png');

%% Torques
%
figure;
hold on;
plot(t, tau_data(1,:), 'LineWidth', 2, 'DisplayName', '$\tau_1$')
hold on;
plot(t, tau_data(2,:), 'LineWidth', 2, 'DisplayName', '$\tau_2$')
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('Torque [Nm]', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
grid on;
set(gca, 'FontSize', 14);
axis tight;
saveas(gcf, 'images/f.png');

%%  Normalized angular velocity
% Need to plot the angular velocity derived by balance of angular momentum
figure;
hold on;
plot(t, theta_d/w, 'LineWidth', 2, 'DisplayName', '$\frac{\dot{\theta}(t)}{\omega^2}$')
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('Normalized angular velocity $[\frac{1/s}{1/s}]$', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
grid on;
set(gca, 'FontSize', 14);
axis tight;
saveas(gcf, 'images/e.png');

%%  Normalized angular acceleration
%
figure;
hold on;
plot(t, theta_dd/w, 'LineWidth', 2, 'DisplayName', '$\frac{\ddot{\theta}(t)}{\omega^2}$')
xlabel('Time (t) [s]', 'Interpreter', 'latex');
ylabel('Normalized angular acceleration $[\frac{1/s^2}{1/s^2}]$', 'Interpreter', 'latex');
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
grid on;
set(gca, 'FontSize', 14);
axis tight;
saveas(gcf, 'images/g.png');

%% Animation
animation(t,X);