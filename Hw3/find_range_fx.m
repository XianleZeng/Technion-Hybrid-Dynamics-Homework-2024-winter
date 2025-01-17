function [fx_min,fx_max] = find_range_fx(r,alpha,mu)

m = 0.1;
g = 10;
x_c = 1;
y_c = 5;

mu1 = mu(1);
mu2 = mu(2);

alpha1 = alpha(1);
alpha2 = alpha(2);

n1x = cos(alpha1);
n1y = sin(alpha1);
n2x = cos(alpha2);
n2y = sin(alpha2);
t1x = cos(alpha1- pi/2);
t1y = sin(alpha1- pi/2);
t2x = cos(alpha2- pi/2);
t2y = sin(alpha2- pi/2);

r1x = r(1,1);
r1y = r(2,1);
r2x = r(1,2);
r2y = r(2,2);

A = [
    t1x - mu1*n1x, t1y - mu1*n1y, 0, 0, 0;
   -t1x - mu1*n1x, -t1y - mu1*n1y, 0, 0, 0;
    0, 0, t2x - mu2*n2x, t2y - mu2*n2y, 0;
    0, 0, -t2x - mu2*n2x, -t2y - mu2*n2y, 0
];

b = [0; 0; 0; 0];

A_tilde = [
    1, 0, 1, 0, 1;
    0, 1, 0, 1, 0;
    -r1y, r1x, -r2y, r2x, -y_c;
];

b_tilde = [0; m*g; m*g*x_c];


c = [0; 0; 0; 0; 1]; 

[X_min, ~, exitflag_min] = linprog(c, A, b, A_tilde, b_tilde);
[X_max, ~, exitflag_max] = linprog(-c, A, b, A_tilde, b_tilde);


if exitflag_min == -2 || exitflag_max == -2
    disp("No feasible solution.")
    fx_min = Nan;
    fx_max = Nan;
    return;
else
    if exitflag_min == 1
        fx_min = [0,0,0,0,1]*X_min;
    elseif exitflag_min == -3 
        fx_min = -inf;
    end

    if exitflag_max == 1
        fx_max = [0,0,0,0,1]*X_max;
    elseif exitflag_max == -3 
        fx_max = inf;
    end
end

%%%%%%%%%%% plot result with force directions %%%%%%%%%%%%%%
figure;
legends = [];

% Contact Normal
line([r1x, r1x + 5*cos(alpha1)], [r1y, r1y + 5*sin(alpha1)], 'Color', 'blue', 'LineStyle', '--', 'LineWidth', 2);
hold on;
line([r2x, r2x + 5*cos(alpha2)], [r2y, r2y + 5*sin(alpha2)], 'Color', 'blue', 'LineStyle', '--', 'LineWidth', 2);

legends = [legends, "Contact Normal", ""];

% Contact Points
plot(r1x, r1y, 'k*', 'LineWidth', 3);
plot(r2x, r2y, 'k*', 'LineWidth', 3);

legends = [legends, "Contact Points", ""];

% Ground Lines
plot([r1x - 0.7*t1x, r1x + 0.7*t1x], [r1y - 0.7*t1y, r1y + 0.7*t1y], 'k-', 'LineWidth', 2);
plot([r2x - 0.7*t2x, r2x + 0.7*t2x], [r2y - 0.7*t2y, r2y + 0.7*t2y], 'k-', 'LineWidth', 2);

legends = [legends, "", ""];

% Friction Cone Lines
line([r1x, r1x + 5*cos(alpha1 + atan(mu1))], [r1y, r1y + 5*sin(alpha1 + atan(mu1))], 'Color', 'red', 'LineWidth', 1.5);
line([r1x, r1x + 5*cos(alpha1 - atan(mu1))], [r1y, r1y + 5*sin(alpha1 - atan(mu1))], 'Color', 'red', 'LineWidth', 1.5);
line([r2x, r2x + 5*cos(alpha2 + atan(mu2))], [r2y, r2y + 5*sin(alpha2 + atan(mu2))], 'Color', 'red', 'LineWidth', 1.5);
line([r2x, r2x + 5*cos(alpha2 - atan(mu2))], [r2y, r2y + 5*sin(alpha2 - atan(mu2))], 'Color', 'red', 'LineWidth', 1.5);

legends = [legends, "Boundary lines of the friction cones", "", "", ""];

% Plot the Center of Mass
plot(x_c, y_c, 'go', 'MarkerSize', 8, 'LineWidth', 2);
legends = [legends, "Center of Mass"];

% Plot resultant force direction range (filled region)

if ~isnan(fx_min) && ~isnan(fx_max)
    % Calculate resultant force magnitudes
    F_min = sqrt(fx_min^2 + (m * g)^2);
    F_max = sqrt(fx_max^2 + (m * g)^2);
    
    % Calculate directions for fx_min and fx_max
    theta_min = atan2(-m * g, fx_min);
    theta_max = atan2(-m * g, fx_max);
    
    % Compute the points of the filled area
    x_fill = [x_c, x_c + F_min * cos(theta_min), x_c + F_max * cos(theta_max)];
    y_fill = [y_c, y_c + F_min * sin(theta_min), y_c + F_max * sin(theta_max)];
    
    % Plot the filled region
    fill(x_fill, y_fill, 'magenta', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on;
    
    % Plot the boundary arrows
    quiver(x_c, y_c, F_min * cos(theta_min), F_min * sin(theta_min), ...
           'Color', 'magenta', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
    quiver(x_c, y_c, F_max * cos(theta_max), F_max * sin(theta_max), ...
           'Color', 'magenta', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
    
    legends = [legends, "Resultant force range"];
end

axis equal;
xlabel('${x}$ [m]', 'Interpreter', 'Latex');
ylabel('${y}$ [m]', 'Interpreter', 'Latex');

legend(legends, 'Interpreter', 'Latex');
title("\textbf{Friction coefficients: $\mu_1=" + num2str(mu1) + ",  \mu_2=" + num2str(mu2) + "$} \ " + ...
      "\textbf{Permissiable force: $" + num2str(fx_min) + " \leq f_x \leq " + num2str(fx_max) + "$}", ...
      'Interpreter', 'latex');
saveas(gcf, ".\images\find_max_fx.png");
end