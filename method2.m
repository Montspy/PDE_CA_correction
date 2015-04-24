
clear

tf = 1;
dt = 0;
dx = 1;

I = imread('test images/axial_numbers.jpg');
R = double(I(:,:,1))/255;
G = double(I(:,:,2))/255;
B = double(I(:,:,3))/255;

alpha = 200;
beta = 1;

S = size(R);
N = S(1);
M = S(2);

figure(1);
image(double(I)/255);
title('Original image');

    % x and y offset
u = zeros(size(R));
v = zeros(size(R));

b = zeros(size(R)); % b = -(R - G)*R_x
c = zeros(size(R)); % c = -(R - G)*R_y

E = [];
tt = [];

tic;

%% Red
t = 0;
while(t < tf)
    % Compute b(x,y) and c(x,y)
    for x = 1:M
        for y = 1:N
            x_off = max(1, min(M, x + round(u(y,x))));
            y_off = max(1, min(N, y + round(v(y,x))));
            
            % Limit indexes
            x_off_p1 = max(1, min(M, x+1 + round(u(y,x))));
            x_off_m1 = max(1, min(M, x-1 + round(u(y,x))));
            y_off_p1 = max(1, min(N, y+1 + round(v(y,x))));
            y_off_m1 = max(1, min(N, y-1 + round(v(y,x))));
            
            x_p1 = min(M, x+1);
            x_m1 = max(1, x-1);
            y_p1 = min(N, y+1);
            y_m1 = max(1, y-1);
            
            % Spatial derivatives (centered)
            R_x = (R(y_off, x_off_p1) - R(y_off, x_off_m1))/(2*dx);

            R_y = (R(y_off_p1, x_off) - R(y_off_m1, x_off))/(2*dx);
            
            % b, c
            b(y,x) = -(R(y_off,x_off) - G(y,x))*R_x;
            c(y,x) = -(R(y_off,x_off) - G(y,x))*R_y;
        end
    end
    
    % Compute maximum dt
    dt = dx^2/4;            % Heat equation for u
    dt = min(dt, dx^2/4);   % Heat equation for v
    
    dt = dt/2;
    
    new_u = zeros(size(u));
    new_v = zeros(size(v));
    % Evolve u and v
    for x = 1:M
        for y = 1:N
            x_p1 = min(M, x+1);
            x_m1 = max(1, x-1);
            y_p1 = min(N, y+1);
            y_m1 = max(1, y-1);
            
            u_xx = (u(y, x_p1) - 2*u(y,x) + u(y,x_m1))/(dx^2);
            u_yy = (u(y_p1, x) - 2*u(y,x) + u(y_m1,x))/(dx^2);
            
            new_u(y,x) = u(y,x) + dt*(alpha*b(y,x) + beta*(u_xx + u_yy));
            
            v_xx = (v(y, x_p1) - 2*v(y,x) + v(y,x_m1))/(dx^2);
            v_yy = (v(y_p1, x) - 2*v(y,x) + v(y_m1,x))/(dx^2);
            
            new_v(y,x) = v(y,x) + dt*(alpha*c(y,x) + beta*(v_xx + v_yy));
        end
    end
    u = new_u;
    v = new_v;
    
    E = [E, energy2(R, G, u, v, dx, alpha, beta)];
    t = t + dt;
    tt = [tt t];
end

J = zeros(N, M, 3);
R_off = zeros(size(R));
for x = 1:M
    for y = 1:N
        x_off = max(1, min(M, x + round(u(y,x))));
        y_off = max(1, min(N, y + round(v(y,x))));
        
        R_off(y,x) = R(y_off, x_off);
    end
end

%% Blue
t = 0;

while(t < tf)
    % Compute b(x,y) and c(x,y)
    for x = 1:M
        for y = 1:N
            x_off = max(1, min(M, x + round(u(y,x))));
            y_off = max(1, min(N, y + round(v(y,x))));
            
            % Limit indexes
            x_off_p1 = max(1, min(M, x+1 + round(u(y,x))));
            x_off_m1 = max(1, min(M, x-1 + round(u(y,x))));
            y_off_p1 = max(1, min(N, y+1 + round(v(y,x))));
            y_off_m1 = max(1, min(N, y-1 + round(v(y,x))));
            
            x_p1 = min(M, x+1);
            x_m1 = max(1, x-1);
            y_p1 = min(N, y+1);
            y_m1 = max(1, y-1);
            
            % Spatial derivatives (centered)
            B_x = (B(y_off, x_off_p1) - B(y_off, x_off_m1))/(2*dx);

            B_y = (B(y_off_p1, x_off) - B(y_off_m1, x_off))/(2*dx);
            
            % b, c
            b(y,x) = -(B(y_off,x_off) - B(y,x))*B_x/255;
            c(y,x) = -(B(y_off,x_off) - B(y,x))*B_y/255;
        end
    end
    
    % Compute maximum dt
    dt = dx^2/4;            % Heat equation for u
    dt = min(dt, dx^2/4);   % Heat equation for v
    
    dt = dt/2;
    
    new_u = zeros(size(u));
    new_v = zeros(size(v));
    % Evolve u and v
    for x = 1:M
        for y = 1:N
            x_p1 = min(M, x+1);
            x_m1 = max(1, x-1);
            y_p1 = min(N, y+1);
            y_m1 = max(1, y-1);
            
            u_xx = (u(y, x_p1) - 2*u(y,x) + u(y,x_m1))/(dx^2);
            u_yy = (u(y_p1, x) - 2*u(y,x) + u(y_m1,x))/(dx^2);
            
            new_u(y,x) = u(y,x) + dt*(b(y,x) + (u_xx + u_yy));
            
            v_xx = (v(y, x_p1) - 2*v(y,x) + v(y,x_m1))/(dx^2);
            v_yy = (v(y_p1, x) - 2*v(y,x) + v(y_m1,x))/(dx^2);
            
            new_v(y,x) = v(y,x) + dt*(c(y,x) + (v_xx + v_yy));
        end
    end
    u = new_u;
    v = new_v;
    
    E = [E, energy2(B, G, u, v, dx, alpha, beta)];
    t = t + dt;
    tt = [tt t];
end

B_off = zeros(size(B));
for x = 1:M
    for y = 1:N
        x_off = max(1, min(M, x + round(u(y,x))));
        y_off = max(1, min(N, y + round(v(y,x))));
        
        B_off(y,x) = B(y_off, x_off);
    end
end

J(:,:,1) = R_off;
J(:,:,2) = G;
J(:,:,3) = B_off;   % TODO

toc

figure(2);
image(J);
title('Corrected image');
