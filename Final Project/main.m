clc
close all
clear

Z_guess = [0; 0; -0.166539; -0.166539; 0; 0; 0.803329; -0.448996]; % Initial guess from reference thesis
Z_star = find_fixed_point(Z_guess);
% Z_star = [0; 0; -0.166539; -0.166539; 0; 0; 0.803329; -0.448996];
q0 = Z_star;
t0 = 0;
time_span = 10; % sec
Ie_stick = 1;
t_stick_data = [];
q_stick_data = [];


while ~isempty(Ie_stick)
op_stick = odeset('reltol',1e-8,'abstol',1e-8,'Events',@(t, q)events_stick(t, q));
[t_stick, q_stick, Te_stick, Xe_stick, Ie_stick] = ode45(@(t,q)sys_stick(t, q), [t0, time_span], q0, op_stick);
t_stick_data = [t_stick_data; t_stick];
q_stick_data = [q_stick_data; q_stick];
    if ~isempty(Xe_stick)
        q0 = impact_law(Xe_stick(end,:).');
        t0 = Te_stick(end);
    end
end

figure;
plot(t_stick_data, q_stick_data(:,3), 'b', 'LineWidth', 1.5); hold on;
plot(t_stick_data, q_stick_data(:,4), 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angles \theta_1, \theta_2 (rad)');
legend('\theta_1(t)', '\theta_2(t)');
title('Angles vs Time');
grid on;
% scatter(Te_stick(2), Xe_stick(2,3), 'k', 'filled'); % 假设 Xe_stick(:,1) 代表 scuffing 事件
