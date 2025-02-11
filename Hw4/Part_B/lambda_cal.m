function [lam_n, lam_t] = lambda_cal(t, q, q_d, t_contact_type, contact_type, separation)

    [m_1, m_2, l, h, J_1, J_2, R, g, mu]=model_params();
    lam_n = zeros(1,length(t));
    lam_t = zeros(1,length(t));
    
    ch = 1;
    for i = 1:1:length(t)
        
        if (length(t_contact_type) >= ch & t_contact_type(ch) <= t(i)) 
            if i == length(t)
                if separation == 1
                    state = 4;
                else
                    state = contact_type(ch);
                end
            else
                state = contact_type(ch);
                ch = ch + 1;
            end
        end
        
                
        if (state == 1)

            [~, lambda] = dyn_sol_stick(t,q(:,i),q_d(:,i));  % A*[q_dd; lambda] = B
            
            lam_t(i) = lambda(1);
            lam_n(i) = lambda(2);
                
        elseif (state == 2)
            slip_dir = 1;
            sigma = 1;
            [~, lambda]=dyn_sol_slip(t,q(:,i),q_d(:,i),slip_dir);
            lam_n(i) = lambda;
            lam_t(i) = -sigma*mu*lam_n(i);
            
        elseif (state == 3)
            slip_dir = 2;
            sigma = -1;
            [~, lambda]=dyn_sol_slip(t,q(:,i),q_d(:,i),slip_dir);
            lam_n(i) = lambda;
            lam_t(i) = -sigma*mu*lam_n(i);
            
        elseif (state == 4)
            lam_n(i) = 0;
            lam_t(i) = 0;
        end
    end
end


