clear

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
tspan = 0:0.01:40;  
X0 = [0; 0; 0; 0; 0; 0];             
[t, X] = ode45(@state_eq, tspan, X0, options); 

x = X(:, 1);
y = X(:, 2);
theta = X(:, 3);

phi_data = [];
for i = 1:length(t)
    [phi, phi_d, phi_dd] = angles_input(t(i)); 
    phi_data = [phi_data, phi];
end
phi1 = phi_data(1,:);
phi2 = phi_data(2,:);

plot(t, x);
hold on;
plot(t, theta);
% hold on;
% plot(t, phi1);
% hold on;
% plot(t, phi2);

xlabel('Time (t)');
ylabel('y(t)');
title('Solution');
grid on;


%%
L = 0.1;  

% animation
figure;
hold on;
grid on;

% 初始化绘图对象
h_center = plot(x(1), y(1), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % 中间杆中心点
h_rod0 = plot([0, 0], [0, 0], 'b-', 'LineWidth', 2); % 中间杆
h_rod1 = plot([0, 0], [0, 0], 'g-', 'LineWidth', 2); % 左杆
h_rod2 = plot([0, 0], [0, 0], 'r-', 'LineWidth', 2); % 右杆


% 动画循环
for i = 1:length(t)
    % 中间杆两端点坐标
    x0_start = x(i) - (L / 2) * cos(theta(i));
    y0_start = y(i) - (L / 2) * sin(theta(i));
    x0_end = x(i) + (L / 2) * cos(theta(i));
    y0_end = y(i) + (L / 2) * sin(theta(i));

    % 左杆端点坐标
    x1 = x0_start - L * cos(theta(i) + phi1(i));
    y1 = y0_start - L * sin(theta(i) + phi1(i));

    % 右杆端点坐标
    x2 = x0_end + L * cos(theta(i) + phi2(i));
    y2 = y0_end + L * sin(theta(i) + phi2(i));

    % 更新图像
    set(h_center, 'XData', x(i), 'YData', y(i));
    set(h_rod0, 'XData', [x0_start, x0_end], 'YData', [y0_start, y0_end]);
    set(h_rod1, 'XData', [x0_start, x1], 'YData', [y0_start, y1]);
    set(h_rod2, 'XData', [x0_end, x2], 'YData', [y0_end, y2]);

    % 动态调整图形范围
    all_x = [x0_start, x0_end, x1, x2];
    all_y = [y0_start, y0_end, y1, y2];
    
    % 计算当前所有点的最小和最大坐标
    x_min = min(all_x) - 0.1;
    x_max = max(all_x) + 0.1;
    y_min = min(all_y) - 0.1;
    y_max = max(all_y) + 0.1;
    
    % 设置动态范围
    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    axis equal
    % 动态刷新
    drawnow;
    pause(0.0005); % 控制动画速度
end