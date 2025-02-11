function [m_1, m_2, l, h, J_1, J_2, R, g, mu] = model_params()

m_1 = 5;
m_2 = 15;
l = 0.2;
h = 0.4;
R = 0.6;
J_1 = 0.1*m_1*R^2;
J_2 = 0.5*m_2*l^2;
g = 10;
mu = 0.05;

end

