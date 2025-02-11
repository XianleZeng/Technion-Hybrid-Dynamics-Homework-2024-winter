function [x_min,x_max] = fric_eq(r,alpha,mu)

m = 0.1;
g = 10;

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
    1, 0, 1, 0, 0;
    0, 1, 0, 1, 0;
    -r1y, r1x, -r2y, r2x, -m*g;
];

b_tilde = [0; m*g; 0];


c = [0; 0; 0; 0; 1]; 

options = optimoptions('linprog', ...
    'OptimalityTolerance', 1e-12, ... % 设置最优性容差
    'ConstraintTolerance', 1e-12, ... % 设置约束容差
    'Display', 'iter'); % 显示迭代过程


[X_min, ~, exitflag_min] = linprog(c, A, b, A_tilde, b_tilde, [], [], options);
[X_max, ~, exitflag_max] = linprog(-c, A, b, A_tilde, b_tilde, [], [], options);


if exitflag_min == -2 || exitflag_max == -2
    disp("No feasible solution.")
    x_min = Nan;
    x_max = Nan;
    return;
else
    if exitflag_min == 1
        x_min = [0,0,0,0,1]*X_min;
    elseif exitflag_min == -3 
        x_min = -inf;
    end

    if exitflag_max == 1
        x_max = [0,0,0,0,1]*X_max;
    elseif exitflag_max == -3 
        x_max = inf;
    end
end


%%%%%%%%%%% plot reslut %%%%%%%%%%%%%%

legends =[];

% Contact Normal
line([r1x,r1x+5*cos(alpha1)],[r1y,r1y+5*sin(alpha1)], 'Color', 'blue', 'LineStyle', '--', 'LineWidth', 2);
hold on;
line([r2x,r2x+5*cos(alpha2)],[r2y,r2y+5*sin(alpha2)], 'Color', 'blue', 'LineStyle', '--', 'LineWidth', 2);
 
legends = [legends, "Contact Normal", ""];
 
% Contact Points
plot(r1x,r1y,'k*','LineWidth',3);
plot(r2x,r2y,'k*','LineWidth',3);
 
legends = [legends, "Contact Points", ""];
 
% Ground Lines
plot([r1x-0.7*t1x, r1x+0.7*t1x], [r1y-0.7*t1y, r1y+0.7*t1y], 'k-', 'LineWidth', 2);
plot([r2x-0.7*t2x, r2x+0.7*t2x], [r2y-0.7*t2y, r2y+0.7*t2y], 'k-', 'LineWidth', 2);
legends = [legends, "", ""];
 
if exitflag_min ~= -2 && exitflag_max ~= -2

    % Feasible Region:
    y_limit = ylim;
    x_limit = xlim;

    x_min_plot = x_min;
    x_max_plot = x_max;
    
    if exitflag_min == -3 % unbounded
        x_min_plot = x_limit(1)-3;
        xlim(x_limit);
    else
        line([x_min_plot,x_min_plot], [y_limit(1)-2,y_limit(2)+2], 'Color', '#9932CC', 'LineWidth', 3);
        legends = [legends, ""];
    end

    if exitflag_max == -3 % unbounded
        x_max_plot = x_limit(2)+10;
        xlim(x_limit);
    else
        line([x_max_plot,x_max_plot], [y_limit(1)-2,y_limit(2)+2], 'Color', '#9932CC', 'LineWidth', 3);
        legends = [legends, ""];
    end    

    fill([x_min_plot,x_min_plot,x_max_plot,x_max_plot,x_min_plot], [y_limit(1)-2,y_limit(2)+2,y_limit(2)+2,y_limit(1)-2,y_limit(1)-2], 'cyan', 'FaceAlpha', 0.3, 'LineStyle', 'none')
    legends = [legends, "﻿COM equilibrium region"];
end

% Friction Cone Lines
line([r1x,r1x+5*cos(alpha1+atan(mu1))],[r1y,r1y+5*sin(alpha1+atan(mu1))], 'Color', 'red', 'LineWidth', 1.5);
line([r1x,r1x+5*cos(alpha1-atan(mu1))],[r1y,r1y+5*sin(alpha1-atan(mu1))], 'Color', 'red', 'LineWidth', 1.5);
line([r2x,r2x+5*cos(alpha2+atan(mu2))],[r2y,r2y+5*sin(alpha2+atan(mu2))], 'Color', 'red', 'LineWidth', 1.5);
line([r2x,r2x+5*cos(alpha2-atan(mu2))],[r2y,r2y+5*sin(alpha2-atan(mu2))], 'Color', 'red', 'LineWidth', 1.5);
 
legends = [legends, "﻿Boundary lines of the friction cones", "", "", ""];


axis equal
xlabel('${x}$ [m]', 'Interpreter', 'Latex')
ylabel('${y}$ [m]', 'Interpreter', 'Latex')
 
legend(legends)
title("\textbf{Friction coefficients: $\mu_1=" + num2str(mu1) + ",  \mu_2=" + num2str(mu2) + "$} \ " + ...
      "\textbf{Permission COM region: $" + num2str(x_min) + " \leq x_c \leq " + num2str(x_max) + "$}", ...
      'Interpreter', 'latex');
saveas(gcf, ".\images\mu1_"+ num2str(mu1) + "mu2_" + num2str(mu2) + ".png");
end

