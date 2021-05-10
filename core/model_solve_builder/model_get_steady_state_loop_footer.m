                y = g(x);
                % increments
                incr_x = norm(x-x_old)/max(1e-10, norm(x));
                incr_y = norm((y-y_old)./y_norm);
                if incr_x/dt < opt.eps_x && incr_y/dt < opt.eps_y
                    steady_state = 1;
                end
            end