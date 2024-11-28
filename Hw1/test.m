t = 0:0.01:10;  
phi_data = [];
for i = 1:length(t)
    [phi, phi_d, phi_dd] = angles_input(t(i)); 
    phi_data = [phi_data, phi];
end
phi1 = phi_data(1,:);
phi2 = phi_data(2,:);

plot(t, phi1);
hold on;
plot(t, phi2);

xlabel('Time (t)');
ylabel('phi');
grid on;

