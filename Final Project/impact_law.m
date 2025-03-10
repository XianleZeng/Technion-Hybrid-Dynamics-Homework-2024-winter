function Xnew = impact_law(Xold)

    [m, m_h, l, I_c, g, alpha, mu]=model_params();
    % Extract pre-impact state
    x = Xold(1);
    y = Xold(2);
    theta_1 = Xold(3);
    theta_2 = Xold(4);
    dx = Xold(5);
    dy = Xold(6);
    dtheta_1 = Xold(7);
    dtheta_2 = Xold(8);
    
    % % Compute post-impact velocities using Chatterjee’s impact law
    % Relabel coordinates
    Xnew = [x + 2*l*sin(theta_1) + 2*l*sin(theta_2);
               y + 2*l*cos(theta_1) - 2*l*cos(theta_2); 
               -theta_2; 
               -theta_1; 
               0; 
               0; 
               -dtheta_2; 
               -dtheta_1];
    
    % % Check for double-foot impact
    % if Xnew(3) < 0
    %     Xnew = 'failure';
    % end
end