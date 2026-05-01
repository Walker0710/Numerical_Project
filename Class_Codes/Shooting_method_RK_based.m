
clc; clear; close all;

% PROBLEM DEFINITION
% y'' = 6y^2 - x

f  = @(x,y,z) 6*y.^2 - x;
fy = @(x,y,z) 12*y;
fz = @(x,y,z) 0;

a = 0; b = 1;
alpha = 1; beta = 5;

h = 1/3;
x = a:h:b;

% NEWTON SHOOTING METHOD

s = 1.2;              % initial guess for y'(0)
tol = 1e-8;

for iter = 1:10

    % Solve IVP using RK4
    [y, eta] = RK4_solver(s, x, h, alpha, f, fy, fz);

    % Residual
    phi  = y(end) - beta;
    dphi = eta(end);

   %%%Safety check 
    if abs(dphi) < 1e-12
       error('Newton method failed: derivative too small');
    end

    % Newton update
    s_new = s - phi/dphi;

    if abs(s_new - s) < tol
        break;
    end

    s = s_new;
end

% RESULTS

fprintf('Final slope s = %.6f\n', s);
fprintf('y(1/3) = %.6f\n', y(2));
fprintf('y(2/3) = %.6f\n', y(3));

disp('   x        y(x)')
disp([x' y'])

% RK4 FUNCTION

function [y, eta] = RK4_solver(s, x, h, alpha, f, fy, fz)

    N = length(x);

    y1 = zeros(1,N); y2 = zeros(1,N);
    e1 = zeros(1,N); e2 = zeros(1,N);

    % Initial conditions
    y1(1) = alpha;
    y2(1) = s;

    e1(1) = 0;
    e2(1) = 1;

    for i = 1:N-1
        xi = x(i);

        % -------- k1 --------
        k1y1 = y2(i);
        k1y2 = f(xi,y1(i),y2(i));

        k1e1 = e2(i);
        k1e2 = fy(xi,y1(i),y2(i))*e1(i) + fz(xi,y1(i),y2(i))*e2(i);

        % -------- k2 --------
        y1_2 = y1(i) + h*k1y1/2;
        y2_2 = y2(i) + h*k1y2/2;
        e1_2 = e1(i) + h*k1e1/2;
        e2_2 = e2(i) + h*k1e2/2;

        k2y1 = y2_2;
        k2y2 = f(xi+h/2, y1_2, y2_2);

        k2e1 = e2_2;
        k2e2 = fy(xi+h/2,y1_2,y2_2)*e1_2 + fz(xi+h/2,y1_2,y2_2)*e2_2;

        % -------- k3 --------
        y1_3 = y1(i) + h*k2y1/2;
        y2_3 = y2(i) + h*k2y2/2;
        e1_3 = e1(i) + h*k2e1/2;
        e2_3 = e2(i) + h*k2e2/2;

        k3y1 = y2_3;
        k3y2 = f(xi+h/2, y1_3, y2_3);

        k3e1 = e2_3;
        k3e2 = fy(xi+h/2,y1_3,y2_3)*e1_3 + fz(xi+h/2,y1_3,y2_3)*e2_3;

        % -------- k4 --------
        y1_4 = y1(i) + h*k3y1;
        y2_4 = y2(i) + h*k3y2;
        e1_4 = e1(i) + h*k3e1;
        e2_4 = e2(i) + h*k3e2;

        k4y1 = y2_4;
        k4y2 = f(xi+h, y1_4, y2_4);

        k4e1 = e2_4;
        k4e2 = fy(xi+h,y1_4,y2_4)*e1_4 + fz(xi+h,y1_4,y2_4)*e2_4;

        % -------- UPDATE --------
        y1(i+1) = y1(i) + h/6*(k1y1 + 2*k2y1 + 2*k3y1 + k4y1);
        y2(i+1) = y2(i) + h/6*(k1y2 + 2*k2y2 + 2*k3y2 + k4y2);

        e1(i+1) = e1(i) + h/6*(k1e1 + 2*k2e1 + 2*k3e1 + k4e1);
        e2(i+1) = e2(i) + h/6*(k1e2 + 2*k2e2 + 2*k3e2 + k4e2);
    end

    y = y1;
    eta = e1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%using interpolation
% clc; clear; close all;
% 
% %% =========================================
% % PROBLEM DEFINITION
% % y'' = 6y^2 - x
% %% =========================================
% 
% f  = @(x,y,z) 6*y.^2 - x;
% fy = @(x,y,z) 12*y;
% fz = @(x,y,z) 0;
% 
% a = 0; b = 1;
% alpha = 1; beta = 5;
% 
% h = 0.01;                 % step size
% x = a:h:b;
% 
% %% =========================================
% % NEWTON SHOOTING METHOD
% %% =========================================
% 
% s = 1.2;                  % initial guess
% tol = 1e-8;
% max_iter = 20;
% 
% for iter = 1:max_iter
% 
%     % Solve IVP using RK4
%     [y, eta] = RK4_solver(s, x, h, alpha, f, fy, fz);
% 
%     % Residual
%     phi  = y(end) - beta;
%     dphi = eta(end);
% 
%     % Safety check 
%     if abs(dphi) < 1e-12
%         error('Newton method failed: derivative too small');
%     end
% 
%     % Newton Raphson update
%     s_new = s - phi/dphi;
% 
%     if abs(s_new - s) < tol
%         s = s_new;
%         break;
%     end
% 
%     s = s_new;
% end
% 
% %% =========================================
% % RESULTS
% %% =========================================
% 
% % Correct evaluation using interpolation
% y_13 = interp1(x, y, 1/3);
% y_23 = interp1(x, y, 2/3);
% 
% fprintf('Final slope s = %.6f\n', s);
% fprintf('y(1/3) = %.6f\n', y_13);
% fprintf('y(2/3) = %.6f\n', y_23);
% 
% disp('   x        y(x)')
% disp([x' y'])
% 
% %% =========================================
% % RK4 FUNCTION
% %% =========================================
% 
% function [y, eta] = RK4_solver(s, x, h, alpha, f, fy, fz)
% 
% N = length(x);
% 
% y = zeros(1,N); z = zeros(1,N);
% e1 = zeros(1,N); e2 = zeros(1,N);
% 
% % Initial conditions
% y(1) = alpha;
% z(1) = s;
% 
% e1(1) = 0;
% e2(1) = 1;
% 
% for i = 1:N-1
%     xi = x(i);
% 
%     % ---- k1 ----
%     k1y = z(i);
%     k1z = f(xi, y(i), z(i));
% 
%     k1e1 = e2(i);
%     k1e2 = fy(xi,y(i),z(i))*e1(i) + fz(xi,y(i),z(i))*e2(i);
% 
%     % ---- k2 ----
%     k2y = z(i) + 0.5*h*k1z;
%     k2z = f(xi+h/2, y(i)+0.5*h*k1y, z(i)+0.5*h*k1z);
% 
%     k2e1 = e2(i) + 0.5*h*k1e2;
%     k2e2 = fy(xi+h/2, y(i)+0.5*h*k1y, z(i)+0.5*h*k1z)*(e1(i)+0.5*h*k1e1) ...
%          + fz(xi+h/2, y(i)+0.5*h*k1y, z(i)+0.5*h*k1z)*(e2(i)+0.5*h*k1e2);
% 
%     % ---- k3 ----
%     k3y = z(i) + 0.5*h*k2z;
%     k3z = f(xi+h/2, y(i)+0.5*h*k2y, z(i)+0.5*h*k2z);
% 
%     k3e1 = e2(i) + 0.5*h*k2e2;
%     k3e2 = fy(xi+h/2, y(i)+0.5*h*k2y, z(i)+0.5*h*k2z)*(e1(i)+0.5*h*k2e1) ...
%          + fz(xi+h/2, y(i)+0.5*h*k2y, z(i)+0.5*h*k2z)*(e2(i)+0.5*h*k2e2);
% 
%     % ---- k4 ----
%     k4y = z(i) + h*k3z;
%     k4z = f(xi+h, y(i)+h*k3y, z(i)+h*k3z);
% 
%     k4e1 = e2(i) + h*k3e2;
%     k4e2 = fy(xi+h, y(i)+h*k3y, z(i)+h*k3z)*(e1(i)+h*k3e1) ...
%          + fz(xi+h, y(i)+h*k3y, z(i)+h*k3z)*(e2(i)+h*k3e2);
% 
%     % ---- update ----
%     y(i+1) = y(i) + h/6*(k1y + 2*k2y + 2*k3y + k4y);
%     z(i+1) = z(i) + h/6*(k1z + 2*k2z + 2*k3z + k4z);
% 
%     e1(i+1) = e1(i) + h/6*(k1e1 + 2*k2e1 + 2*k3e1 + k4e1);
%     e2(i+1) = e2(i) + h/6*(k1e2 + 2*k2e2 + 2*k3e2 + k4e2);
% end
% 
% eta = e1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Not useful always
% clc; clear; close all;
% 
% %% =========================================
% % PROBLEM DEFINITION
% % y'' = 6y^2 - x
% %% =========================================
% 
% f = @(x,y,z) 6*y.^2 - x;
% 
% a = 0; b = 1;
% alpha = 1; beta = 5;
% 
% h = 1/3;
% x = a:h:b;
% 
% %% =========================================
% % NEWTON SHOOTING (FINITE DIFFERENCE)
% %% =========================================
% 
% s = 1.2;              % initial guess
% tol = 1e-8;
% 
% for iter = 1:10
% 
%     % Solve IVP
%     y = RK4_solver(s, x, h, alpha, f);
%     F = y(end) - beta;
% 
%     % Finite difference derivative
%     ds = 1e-6;
%     y_pert = RK4_solver(s + ds, x, h, alpha, f);
%     Fp = (y_pert(end) - beta - F) / ds;
% 
%     % Newton update
%     s_new = s - F/Fp;
% 
%     if abs(s_new - s) < tol
%         break;
%     end
% 
%     s = s_new;
% end
% 
% %% =========================================
% % RESULTS
% %% =========================================
% 
% fprintf('Final slope s = %.6f\n', s);
% fprintf('y(1/3) = %.6f\n', y(2));
% fprintf('y(2/3) = %.6f\n', y(3));
% 
% disp('   x        y(x)')
% disp([x' y'])
% 
% %% =========================================
% % RK4 FUNCTION (NO SENSITIVITY)
% %% =========================================
% 
% function y = RK4_solver(s, x, h, alpha, f)
% 
%     N = length(x);
% 
%     y1 = zeros(1,N);
%     y2 = zeros(1,N);
% 
%     % Initial conditions
%     y1(1) = alpha;
%     y2(1) = s;
% 
%     for i = 1:N-1
%         xi = x(i);
% 
%         % -------- k1 --------
%         k1y1 = y2(i);
%         k1y2 = f(xi,y1(i),y2(i));
% 
%         % -------- k2 --------
%         y1_2 = y1(i) + h*k1y1/2;
%         y2_2 = y2(i) + h*k1y2/2;
% 
%         k2y1 = y2_2;
%         k2y2 = f(xi+h/2, y1_2, y2_2);
% 
%         % -------- k3 --------
%         y1_3 = y1(i) + h*k2y1/2;
%         y2_3 = y2(i) + h*k2y2/2;
% 
%         k3y1 = y2_3;
%         k3y2 = f(xi+h/2, y1_3, y2_3);
% 
%         % -------- k4 --------
%         y1_4 = y1(i) + h*k3y1;
%         y2_4 = y2(i) + h*k3y2;
% 
%         k4y1 = y2_4;
%         k4y2 = f(xi+h, y1_4, y2_4);
% 
%         % -------- UPDATE --------
%         y1(i+1) = y1(i) + h/6*(k1y1 + 2*k2y1 + 2*k3y1 + k4y1);
%         y2(i+1) = y2(i) + h/6*(k1y2 + 2*k2y2 + 2*k3y2 + k4y2);
%     end
% 
%     y = y1;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%Bvp-4c solver%%%%%%
% %clc; clear; close all;
% 
% %% INITIAL MESH
% xmesh = linspace(0,1,10);
% 
% % Initial guess (important!)
% solinit = bvpinit(xmesh, @guess);
% 
% %% SOLVE USING bvp4c
% sol = bvp4c(@odefun, @bcfun, solinit);
% 
% %% EXTRACT SOLUTION
% x = sol.x;
% y = sol.y(1,:);
% 
% % Required values
% y_13 = deval(sol, 1/3);
% y_23 = deval(sol, 2/3);
% 
% fprintf('y(1/3) = %.6f\n', y_13(1));
% fprintf('y(2/3) = %.6f\n', y_23(1));
% 
% disp('   x        y(x)')
% disp([x' y'])
% 
% %% =========================================
% % ODE SYSTEM
% %% =========================================
% function dydx = odefun(x,y)
% 
%     dydx = zeros(2,1);
%     dydx(1) = y(2);
%     dydx(2) = 6*y(1)^2 - x;
% end
% 
% %% =========================================
% % BOUNDARY CONDITIONS
% %% =========================================
% function res = bcfun(ya,yb)
% 
%     res = [ ya(1) - 1;    % y(0)=1
%             yb(1) - 5 ];  % y(1)=5
% end
% 
% %% =========================================
% % INITIAL GUESS FUNCTION
% %% =========================================
% function g = guess(x)
% 
%     g = [1 + 4*x;   % simple linear guess
%          4];        % derivative guess
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%Best Approach combined with shooting and BVP-4c solver%%%%%
% clc; clear; close all;
% 
% %% =========================================
% % PROBLEM
% % y'' = 6y^2 - x
% %% =========================================
% 
% f  = @(x,y,z) 6*y.^2 - x;
% fy = @(x,y,z) 12*y;
% fz = @(x,y,z) 0;
% 
% a = 0; b = 1;
% alpha = 1; beta = 5;
% 
% h = 0.01;
% xspan = a:h:b;
% 
% %% =========================================
% % STEP 1: SHOOTING (NEWTON METHOD)
% %% =========================================
% 
% s = 1.5;       % initial guess
% tol = 1e-8;
% 
% for iter = 1:20
% 
%     [y, z, eta] = RK4_NR(s, xspan, h, alpha, f, fy, fz);
% 
%     F  = y(end) - beta;
%     dF = eta(end);   % exact derivative
% 
%     s_new = s - F/dF;
% 
%     if abs(s_new - s) < tol
%         break;
%     end
% 
%     s = s_new;
% end
% 
% fprintf('Newton shooting slope ≈ %.6f\n', s);
% 
% % Store shooting solution
% xs = xspan; ys = y; zs = z;
% 
% %% =========================================
% % STEP 2: bvp4c USING SHOOTING GUESS
% %% =========================================
% 
% solinit = bvpinit(xs, @(x) guess_from_shoot(x, xs, ys, zs));
% 
% sol = bvp4c(@odefun, @bcfun, solinit);
% 
% %% =========================================
% % RESULTS
% %% =========================================
% 
% x = sol.x;
% y = sol.y(1,:);
% 
% y_13 = deval(sol, 1/3);
% y_23 = deval(sol, 2/3);
% 
% fprintf('\nFinal (bvp4c refined):\n');
% fprintf('y(1/3) = %.6f\n', y_13(1));
% fprintf('y(2/3) = %.6f\n', y_23(1));
% 
% disp('   x        y(x)')
% disp([x' y'])
% 
% %% =========================================
% % RK4 + NEWTON
% %% =========================================
% 
% function [y, z, eta] = RK4_NR(s, x, h, alpha, f, fy, fz)
% 
% N = length(x);
% 
% y = zeros(1,N); z = zeros(1,N);
% e1 = zeros(1,N); e2 = zeros(1,N);
% 
% y(1) = alpha;
% z(1) = s;
% 
% e1(1) = 0;
% e2(1) = 1;
% 
% for i = 1:N-1
%     xi = x(i);
% 
%     % ---- k1 ----
%     f1 = z(i);
%     f2 = f(xi,y(i),z(i));
% 
%     g1 = e2(i);
%     g2 = fy(xi,y(i),z(i))*e1(i) + fz(xi,y(i),z(i))*e2(i);
% 
%     % ---- k2 ----
%     f1_2 = z(i) + 0.5*h*f2;
%     f2_2 = f(xi+h/2, y(i)+0.5*h*f1, z(i)+0.5*h*f2);
% 
%     g1_2 = e2(i) + 0.5*h*g2;
%     g2_2 = fy(xi+h/2,y(i)+0.5*h*f1,z(i)+0.5*h*f2)*(e1(i)+0.5*h*g1) ...
%          + fz(xi+h/2,y(i)+0.5*h*f1,z(i)+0.5*h*f2)*(e2(i)+0.5*h*g2);
% 
%     % ---- k3 ----
%     f1_3 = z(i) + 0.5*h*f2_2;
%     f2_3 = f(xi+h/2, y(i)+0.5*h*f1_2, z(i)+0.5*h*f2_2);
% 
%     g1_3 = e2(i) + 0.5*h*g2_2;
%     g2_3 = fy(xi+h/2,y(i)+0.5*h*f1_2,z(i)+0.5*h*f2_2)*(e1(i)+0.5*h*g1_2) ...
%          + fz(xi+h/2,y(i)+0.5*h*f1_2,z(i)+0.5*h*f2_2)*(e2(i)+0.5*h*g2_2);
% 
%     % ---- k4 ----
%     f1_4 = z(i) + h*f2_3;
%     f2_4 = f(xi+h, y(i)+h*f1_3, z(i)+h*f2_3);
% 
%     g1_4 = e2(i) + h*g2_3;
%     g2_4 = fy(xi+h,y(i)+h*f1_3,z(i)+h*f2_3)*(e1(i)+h*g1_3) ...
%          + fz(xi+h,y(i)+h*f1_3,z(i)+h*f2_3)*(e2(i)+h*g2_3);
% 
%     % ---- update ----
%     y(i+1) = y(i) + h/6*(f1 + 2*f1_2 + 2*f1_3 + f1_4);
%     z(i+1) = z(i) + h/6*(f2 + 2*f2_2 + 2*f2_3 + f2_4);
% 
%     e1(i+1) = e1(i) + h/6*(g1 + 2*g1_2 + 2*g1_3 + g1_4);
%     e2(i+1) = e2(i) + h/6*(g2 + 2*g2_2 + 2*g2_3 + g2_4);
% end
% 
% eta = e1;
% end
% 
% %% =========================================
% % GUESS FUNCTION
% %% =========================================
% 
% function g = guess_from_shoot(x, xs, ys, zs)
%     g = [interp1(xs, ys, x);
%          interp1(xs, zs, x)];
% end
% 
% %% =========================================
% % ODE
% %% =========================================
% 
% function dydx = odefun(x,y)
%     dydx = [y(2);
%             6*y(1)^2 - x];
% end
% 
% %% =========================================
% % BC
% %% =========================================
% 
% function res = bcfun(ya,yb)
%     res = [ya(1)-1;
%            yb(1)-5];
% end