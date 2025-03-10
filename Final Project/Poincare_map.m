function Znew = Poincare_map(Zold)
    options = odeset('Events', @events_stick);
    [~, Z, ~, ~, ie] = ode45(@sys_stick, [0 10], Zold, options);
    if ie == 0 
        Znew = 'failure';
    else
        Znew = impact_law(Z(end, :)');
    end
end
