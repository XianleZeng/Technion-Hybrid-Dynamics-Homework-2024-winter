function [phi, phi_d, phi_dd] = input_angle(t)

phi = (pi/2)*cos(t);
phi_d = -(pi/2)*sin(t);
phi_dd = -(pi/2)*cos(t);

end

