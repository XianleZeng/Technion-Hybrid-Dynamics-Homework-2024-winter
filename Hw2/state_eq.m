function X_d=state_eq(t,X)
% STATE_EQ    Dynamic state equation of motion
% Xianle Zeng

%% assign the variables 
qp = X(1:4);  qp_d = X(5:8);

[qb_dd, ~] = dyn_sol(qb,qb_d,t);

X_d = [qb_d; qb_dd];

