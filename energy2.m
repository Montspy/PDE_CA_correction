function [E] = energy2(R, G, u, v, dx, alpha, beta)
    E = 0;
    S = size(R);
    N = S(1);
    M = S(2);

    for x = 1:M
        for y = 1:N
            x_off = max(1, min(M, x + round(u(y,x))));
            y_off = max(1, min(N, y + round(v(y,x))));
            
            % Limit indexes
            x_p1 = min(M, x+1);
            x_m1 = max(1, x-1);
            y_p1 = min(N, y+1);
            y_m1 = max(1, y-1);
            
            % Spatial derivatives (centered)
            u_x = (u(y, x_p1) - u(y, x_m1))/(2*dx);
            u_y = (u(y_p1, x) - u(y_m1, x))/(2*dx);
            
            v_x = (v(y, x_p1) - v(y, x_m1))/(2*dx);
            v_y = (v(y_p1, x) - v(y_m1, x))/(2*dx);
            
            E = E + (alpha*1/2*(R(y_off,x_off) - G(y,x))^2 + beta*1/2*norm([u_x; u_y])^2 + beta*1/2*norm([v_x; v_y])^2)*dx^2;
        end
    end
end

