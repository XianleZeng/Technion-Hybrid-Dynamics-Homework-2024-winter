function animation(t,X)

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

x = X(:, 1);
y = X(:, 2);
theta = X(:, 3);
% theta = rad2deg(theta);
x_d = X(:, 4);
y_d = X(:, 5);
theta_d = X(:, 6);
q = [x, y, theta, phi1, phi2];
q_d = [x_d, y_d, theta_d, phi1_d, phi2_d];

p_cm_data = [];
v_cm_data = [];
for i = 1:length(t)
    [p_cm, v_cm] = center_of_mass(q(i,:), q_d(i,:));
    p_cm_data = [p_cm_data, p_cm];
    v_cm_data = [v_cm_data, v_cm];
end
x_cm = p_cm_data(1,:);
y_cm = p_cm_data(2,:);

L = 0.1;  

% Initialization
figure;
hold on;
grid on;

h_center = plot(x(1), y(1), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % center point of the mid_link
h_CoM = plot(x_cm(1), y_cm(1), 'bo', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'b');
h_rod0 = plot([0, 0], [0, 0], 'b-', 'LineWidth', 2); % mid_link
h_rod1 = plot([0, 0], [0, 0], 'g-', 'LineWidth', 2); % left_link
h_rod2 = plot([0, 0], [0, 0], 'r-', 'LineWidth', 2); % right_link
xlim([-0.2, 0.2]);
for i = 1:length(t)
    % End coordinates of the left link of the middle link
    x0_start = x(i) - (L / 2) * cos(theta(i));
    y0_start = y(i) - (L / 2) * sin(theta(i));
    x0_end = x(i) + (L / 2) * cos(theta(i));
    y0_end = y(i) + (L / 2) * sin(theta(i));

    % End coordinate of the left link
    x1 = x0_start - L * cos(theta(i) + phi1(i));
    y1 = y0_start - L * sin(theta(i) + phi1(i));

    % End coordinate of the right link
    x2 = x0_end + L * cos(theta(i) + phi2(i));
    y2 = y0_end + L * sin(theta(i) + phi2(i));

    % Upate figure
    set(h_center, 'XData', x(i), 'YData', y(i));
    set(h_CoM, 'XData', x_cm(i), 'YData', y_cm(i));
    set(h_rod0, 'XData', [x0_start, x0_end], 'YData', [y0_start, y0_end]);
    set(h_rod1, 'XData', [x0_start, x1], 'YData', [y0_start, y1]);
    set(h_rod2, 'XData', [x0_end, x2], 'YData', [y0_end, y2]);

    all_x = [x0_start, x0_end, x1, x2];
    all_y = [y0_start, y0_end, y1, y2];

    x_min = min(all_x) - 0.1;
    x_max = max(all_x) + 0.1;
    y_min = min(all_y) - 0.1;
    y_max = max(all_y) + 0.1;

    ylim([y(i)-0.2, y(i)+0.2]);

    drawnow;
    pause(0.001); 
end

