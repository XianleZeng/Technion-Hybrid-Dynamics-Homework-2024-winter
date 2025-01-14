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

saveas(gcf, ".\images\mu1_"+ num2str(mu1) + "mu2_" + num2str(mu2) + ".pdf");
end

