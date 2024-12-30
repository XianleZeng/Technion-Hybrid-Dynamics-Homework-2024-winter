function X_d=state_eq_Q4(t,X)
% STATE_EQ    Dynamic state equation of motion
% Xianle Zeng

%% assign the variables 
qp = X(1:3);  qp_d = X(4:6);

[qp_dd, ~] = dyn_sol_Q4(t,qp,qp_d,false);

X_d = [qp_d; qp_dd];
