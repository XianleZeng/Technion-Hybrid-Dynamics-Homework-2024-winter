function Z_star = find_fixed_point(Z_guess)
    P = @(Z) Poincare_map(Z) - Z;
    
    options = optimoptions('fsolve', 'Display', 'iter','FunctionTolerance', 1e-6, 'StepTolerance', 1e-6, 'MaxIterations', 500);
    Z_star = fsolve(P, Z_guess, options);
    
    disp('Fixed Point:');
    disp(Z_star);
end


