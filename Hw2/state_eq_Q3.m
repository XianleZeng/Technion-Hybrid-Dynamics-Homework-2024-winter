function X_d=state_eq_Q3(t,X)
% STATE_EQ    Dynamic state equation of motion
% Xianle Zeng

%% assign the variables 
q = X(1:4);  q_d = X(5:8);

[q_dd, ~] = dyn_sol_Q3(t,q,q_d);  % A*[q_dd; lambda] = B

X_d = [q_d; q_dd]; % X_d = f(X) ; X = [x, y, theta, phi, dx, dy, dtheta, dphi]
