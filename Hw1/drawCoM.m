function drawCoM(center, radius)
    
    theta = linspace(0, 2*pi, 100);
    x = radius * cos(theta) + center(1);
    y = radius * sin(theta) + center(2);
    plot(x, y, 'b', 'LineWidth', 2); hold on;
    
    plot([center(1)-radius, center(1)+radius], [center(2), center(2)], 'r', 'LineWidth', 2); % 水平线
    plot([center(1), center(1)], [center(2)-radius, center(2)+radius], 'r', 'LineWidth', 2); % 垂直线
    
    plot(center(1), center(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(center(1), center(2), ' CoM', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
    
    axis equal;
    grid on;
end