function [phi, phi_d, phi_dd] = input_angle(t)

phi = (pi/4)*cos(t);
phi_d = -(pi/4)*sin(t);
phi_dd = -(pi/4)*cos(t);

end

