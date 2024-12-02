function X_d=state_eq(t,X)
% STATE_EQ    Dynamic state equation of motion
% Input - state vector X = (qb, dqb)
% Output - time-derivative vector X_d(t) = (dqb, ddqb)
% Xianle Zeng
% 27-Nov-2024 13:35:51

%% assign the variables 
qb = X(1:3);  qb_d = X(4:6);

[qb_dd, ~] = dyn_sol(qb,qb_d,t);

X_d = [qb_d; qb_dd];








