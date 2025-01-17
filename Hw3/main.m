clc
close all;
clear;

r1 = [-1; 1];
r2 = [2; 2];
r = [r1,r2];
alpfa = [deg2rad(45); deg2rad(135)];
mu = [0.4; 0.4];
[x_min,x_max] = fric_eq(r,alpfa,mu)
[fx_min,fx_max] = find_range_fx(r,alpfa,mu)

% Define variables
h = 5;
fd_max_data = [];
beta_data = 0:0.01*pi:2*pi;

% Loop to calculate fd_max for each beta
for beta = beta_data
    fd_max = find_max_fd(r, alpfa, mu, beta, h); % Replace with your function definition
    fd_max_data = [fd_max_data, fd_max];
end

% Define variables
h = 5;
fd_max_data = [];
beta_data = 0:0.05:2*pi;

% Loop to calculate fd_max for each beta
for beta = beta_data
    fd_max = find_max_fd(r, alpfa, mu, beta, h); 
    fd_max_data = [fd_max_data, fd_max];
end

% Find the minimum value and its index
[fd_min, idx_min] = min(fd_max_data);
beta_min = beta_data(idx_min);

% Plotting the results
figure;
plot(beta_data, fd_max_data, 'LineWidth', 1.5); % Plot the main curve
hold on;
plot(beta_min, fd_min, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); % Mark the minimum point
hold off;

% Add labels and title
xlabel('\beta (rad)', 'FontSize', 12); % X-axis label
ylabel('Maximum fd', 'FontSize', 12); % Y-axis label
title('h = 5', 'FontSize', 14); % Title
grid on; % Add grid for better visualization
legend('Maximum fd vs \beta', 'Minimum Point', 'FontSize', 12); % Legend

% Annotate the minimum point
text(beta_min, fd_min, ...
    sprintf('  Min: \\beta = %.2f, fd = %.4f', beta_min, fd_min), ...
    'FontSize', 10, 'VerticalAlignment', 'bottom');
saveas(gcf, ".\images\fd_max_h_5.png");



%% Define h_data and initialize variables
h_data = 1:0.1:10;
fd_max_fessible_data = [];
beta_data = 0:0.01*pi:2*pi;

% Loop to calculate the minimum feasible fd_max for each h
for h = h_data
    fd_max_data2 = [];
    for beta = beta_data
        fd_max = find_max_fd(r, alpfa, mu, beta, h); % Replace with your function definition
        fd_max_data2 = [fd_max_data2, fd_max];
    end
    fd_max_fessible = min(fd_max_data2); % Find the minimum fd_max for current h
    fd_max_fessible_data = [fd_max_fessible_data, fd_max_fessible];
end

% Plotting the results
figure;
plot(h_data, fd_max_fessible_data, 'b-', 'LineWidth', 1.5); % Plot data with enhanced style
xlabel('h (Height of COM)', 'FontSize', 12); % X-axis label
ylabel('Minimum Feasible fd', 'FontSize', 12); % Y-axis label
title('Variation of maximum fd with h', 'FontSize', 14); % Title
grid on; % Add grid for better visualization
legend('fd vs h', 'FontSize', 12, 'Location', 'best'); % Legend
hold off;

saveas(gcf, ".\images\fd_max_h_1-10.png");