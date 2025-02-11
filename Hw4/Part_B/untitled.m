function [t_total, q_total, t_contact_type, contact_type, separation] = full_transition_simulation(q0, t0, time_span, contact, slip, separation, L, R, m1, m2, mu, g);
 
    t_total = [];
    q_total = [];
 
    contact_type = [];
    t_contact_type = [0];
    
    Ie_stick = -1;
    Ie_slip = -1;
 
    while (separation == 0 & ~isempty(Ie_stick) & ~isempty(Ie_slip))
 
        if (contact == 1)
            op_stick = odeset('reltol',1e-8,'abstol',1e-8,'Events',@(t, q)events_stick(t, q, L, R, m1, m2, mu, g));
            [t_stick, q_stick, Te_stick, Xe_stick, Ie_stick] = ode45(@(t,q)sys_stick(t, q, L, R, m1, m2, mu, g), [t0, time_span], q0, op_stick);
 
            t_total = [t_total; t_stick];
            q_total = [q_total; q_stick];
 
            t0 = t_stick(end);
            q0 = Xe_stick;
 
            t_contact_type = [t_contact_type; Te_stick];
            contact_type = [contact_type; 1]; % STICK
 
            if ~isempty(Ie_stick)
               contact = 0;
               slip = 1;
               slip_dir = Ie_stick;
            end
        end
 
        if (slip == 1)
 
            op_slip = odeset('reltol',1e-8,'abstol',1e-8,'Events',@(t, q)events_slip(t, q, L, R, m1, m2, mu, g, slip_dir));
            [t_slip, q_slip, Te_slip, Xe_slip, Ie_slip] = ode45(@(t, q)sys_slip(t, q, L, R, m1, m2, mu, g, slip_dir), [t0, time_span], q0, op_slip);
 
            t_total = [t_total; t_slip];
            q_total = [q_total; q_slip];
 
            t0 = t_slip(end);
            q0 = Xe_slip;
 
            t_contact_type = [t_contact_type; Te_slip];
 
            if (slip_dir == 1)
                contact_type = [contact_type; 2]; % SLIP (positive)
            elseif (slip_dir == 2)
                contact_type = [contact_type; 3]; % SLIP (negative)
            end
 
            if (Ie_slip == 1)
 
                th = q_slip(end, 3);
                th_d = q_slip(end, 7);
 
                lam_t =(m1*m2*(4*L*sin(th)*th_d^2 + 3*g*cos(th)*sin(th)))/(2*(6*m1 + 4*m2 - 3*m2*cos(th)^2));
                lam_n = (12*g*m1^2 + 2*g*m2^2 + 11*g*m1*m2 + 2*L*m2^2*th_d^2*cos(th) + 3*g*m1*m2*cos(th)^2 ...
                        + 12*L*m1*m2*th_d^2*cos(th))/(12*m1 + 8*m2 - 6*m2*cos(th)^2);
 
                if (abs(lam_t) < mu*lam_n)
                   contact = 1;
                   slip = 0;
                else
                    if slip_dir == 1
                        slip_dir = 2;
                    elseif slip_dir == 2
                        slip_dir = 1;
                    end                    
                end
 
            elseif (Ie_slip == 2)
                separation = 1;
                t_contact_type = [t_contact_type; Te_slip];
                contact_type = [contact_type; 4]; % STICK
            end
        end
    end
end
