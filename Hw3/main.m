clc
close all;
clear;

r1 = [-1; 1];
r2 = [2; 2];
r = [r1,r2];
alpfa = [deg2rad(45); deg2rad(135)];
mu = [0.4; 0.4];
[x_min,x_max] = fric_eq(r,alpfa,mu)
[fx_min,fx_max] = find_range_fx(r,alpfa,mu)

h = 5;
fd_max_data = 
for beta=0:2*pi
    fd_max = find_max_fd(r,alpfa,mu,beta,h);
    
end