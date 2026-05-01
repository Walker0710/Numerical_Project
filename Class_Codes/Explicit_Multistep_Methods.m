%%%%%Adams Bashforth Methods of third order%%%%%% 
h=0.1;
b=0.5;
m=3;
a=0;
N=int16((b-a)/h);
L = linspace(a,b,N+1);
S = zeros(1,N+1);
% intial condition 
S(1)=1;

% ---- Step 1: Use Taylor series method for first 2 steps ----
for n=1:2                                          
     x=L(n);
     y=S(n);

     y1=x^2+y^2;
     y2=2*x+2*y*y1;
     y3=2+2*(y*y2+y1^2);

     y=y+h*y1+(h^2/2)*y2+(h^3/6)*y3;

     S(n+1)=y;
     fprintf('n= %2.0f, Taylor start value = %12.12f\n\n',n,S(n+1));
end

% ---- Step 2: Adams-Bashforth 3rd Order ----
for j=m:N
S(j+1)=S(j)+(h/12)*(23*(L(j)^2+S(j)^2)-16*(L(j-1)^2+S(j-1)^2)+5*(L(j-2)^2+S(j-2)^2));
fprintf('j= %2.0f, AB3 solution= %12.12f\n\n',j,S(j+1));
end

% ---- Exact solution (symbolic) ----

syms u(t)

usol(t)=dsolve(diff(u,t)==(t^2+u(t)^2),u(0)==1);

% ---- Error ----

Error=abs(usol(b)-S(N+1));
fprintf('Error= %12.12f\n\n',Error);





% %%%%% Adams-Bashforth 3rd Order Method %%%%%
% 
% %clear; clc;
% 
% h = 0.1;
% a = 0;
% b = 0.5;
% 
% N = round((b-a)/h);
% L = linspace(a,b,N+1);
% S = zeros(1,N+1);
% 
% f = @(x,y) x^2 + y^2;
% 
% % Initial condition
% S(1) = 1;
% 
% %% Use RK4 for first 2 steps ----
% for n = 1:2
%     x = L(n);
%     y = S(n);
% 
%     k1 = f(x,y);
%     k2 = f(x + h/2, y + h*k1/2);
%     k3 = f(x + h/2, y + h*k2/2);
%     k4 = f(x + h, y + h*k3);
% 
%     S(n+1) = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
% 
%     fprintf('n= %2d, RK4 start value= %12.12f\n\n', n, S(n+1));
% end
% 
% %% ---- Step 2: Adams-Bashforth 3rd Order ----
% for j = 3:N
%     S(j+1) = S(j) + (h/12)*...
%         (23*f(L(j),S(j)) ...
%         -16*f(L(j-1),S(j-1)) ...
%         +5*f(L(j-2),S(j-2)));
% 
%     fprintf('j= %2d, AB3 Solution= %12.12f\n\n', j, S(j+1));
% end
% 
% %% ---- Exact solution (symbolic) ----
% syms u(t)
% usol(t) = dsolve(diff(u,t) == (t^2 + u(t)^2), u(0) == 1);
% 
% %% Error
% 
% Error = abs((usol(b)) - S(N+1));
% 
% fprintf('Final Error = %12.12f\n\n', Error);







% %%%%% Adams-Bashforth 4th Order Method %%%%%
% 
% clear; clc;
% 
% h = 0.1;
% a = 0;
% b = 0.5;
% 
% N = round((b-a)/h);
% L = linspace(a,b,N+1);
% S = zeros(1,N+1);
% 
% f = @(x,y) x^2 + y^2;
% 
% % Initial condition
% S(1) = 1;
% 
% %% Use RK4 for first 3 steps ----
% for n = 1:3
%     x = L(n);
%     y = S(n);
% 
%     k1 = f(x,y);
%     k2 = f(x + h/2, y + h*k1/2);
%     k3 = f(x + h/2, y + h*k2/2);
%     k4 = f(x + h, y + h*k3);
% 
%     S(n+1) = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
% 
%     fprintf('n= %2d, RK4 start value= %12.12f\n\n', n, S(n+1));
% end
% 
% %% ---- Step 2: Adams-Bashforth 4th Order ----
% for j = 4:N
%     S(j+1) = S(j) + (h/24)*...
%         (55*f(L(j),S(j)) ...
%         -59*f(L(j-1),S(j-1)) ...
%         +37*f(L(j-2),S(j-2)) ...
%         -9*f(L(j-3),S(j-3)));
%     fprintf('j= %2d, AB4 Solution= %12.12f\n\n', j, S(j+1));
% end
% 
% %% ---- Exact solution (symbolic) ----
% syms u(t)
% usol(t) = dsolve(diff(u,t) == (t^2 + u(t)^2), u(0) == 1);
% 
% %% Error
% 
% Error = abs((usol(b)) - S(N+1));
% 
% fprintf('Final Error = %12.12f\n\n', Error);







% %%Predictor Corrector Method (without iteration)
% %%%%% Adams-Bashforth 4 + Adams-Moulton 4 %%%%%
% 
% %clear; clc;
% 
% h = 0.1;
% a = 0;
% b = 0.5;
% 
% N = round((b-a)/h);
% L = linspace(a,b,N+1);
% S = zeros(1,N+1);
% 
% f = @(x,y) x^2 + y^2;
% 
% % Initial condition
% S(1) = 1;
% 
% %% ---- Step 1: RK4 for first 3 steps ----
% for n = 1:3
%     x = L(n);
%     y = S(n);
% 
%     k1 = f(x,y);
%     k2 = f(x + h/2, y + h*k1/2);
%     k3 = f(x + h/2, y + h*k2/2);
%     k4 = f(x + h, y + h*k3);
% 
%     S(n+1) = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
% 
%     fprintf('n= %2d, RK4 start value= %12.12f\n\n', n, S(n+1));
% end
% 
% %% ---- Step 2: AB4 Predictor + AM4 Corrector ----
% for j = 4:N
%     % ---- Predictor (AB4) ----
%     yp = S(j) + (h/24)*...
%         (55*f(L(j),S(j)) ...
%         -59*f(L(j-1),S(j-1)) ...
%         +37*f(L(j-2),S(j-2)) ...
%         -9*f(L(j-3),S(j-3)));
% 
%     % ---- Corrector (AM4) ----
%     S(j+1) = S(j) + (h/24)*...
%         (9*f(L(j+1),yp) ...
%         +19*f(L(j),S(j)) ...
%         -5*f(L(j-1),S(j-1)) ...
%         +f(L(j-2),S(j-2)));
% 
%     fprintf('j= %2d, AB4-AM4 Solution= %12.12f\n\n', j, S(j+1));
% end
% 
% %% ---- Exact solution ----
% syms u(t)
% usol(t) = dsolve(diff(u,t) == (t^2 + u(t)^2), u(0) == 1);
% 
% Error = abs(double(usol(b)) - S(N+1));
% 
% fprintf('Final Error = %12.12f\n\n', Error);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% ---- solution using ode45 ----
% % [x_ref, y_ref] = ode45(f, [a b], 1);
% % y_exact = y_ref(end);
% % 
% % Error = abs(y_exact - S(N+1));
% % 
% % fprintf('Final Error = %12.12e\n\n', Error);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% %%%%%%Predictor Corrector Method %%%%%%%
% 
% %%%%% Adams-Bashforth 4 + Iterated Adams-Moulton 4 %%%%%
% 
% %clear; clc;
% 
% h = 0.1;
% a = 1;
% b = 1.4;
% 
% N = round((b-a)/h);
% L = linspace(a,b,N+1);
% S = zeros(1,N+1);
% 
% f = @(x,y) 1./x.^2 - y./x;
% 
% % Initial condition
% S(1) = 1;
% 
% %% ---- Step 1: RK4 for first 3 steps ----
% for n = 1:3
%     x = L(n);
%     y = S(n);
% 
%     k1 = f(x,y);
%     k2 = f(x + h/2, y + h*k1/2);
%     k3 = f(x + h/2, y + h*k2/2);
%     k4 = f(x + h, y + h*k3);
% 
%     S(n+1) = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
% 
%     fprintf('n= %2d, RK4 start value = %12.12f\n\n', n, S(n+1));
% end
% 
% %% ---- Step 2: AB4 Predictor + Iterated AM4 Corrector ----
% tol = 1e-10;     % convergence tolerance
% max_iter = 10;   % limit
% 
% for j = 4:N
% 
%     % ---- Predictor (AB4) ----
%     yp = S(j) + (h/24)*...
%         (55*f(L(j),S(j)) ...
%         -59*f(L(j-1),S(j-1)) ...
%         +37*f(L(j-2),S(j-2)) ...
%         -9*f(L(j-3),S(j-3)));
% 
%     % ---- Iterative Corrector (AM4) ----
%     yc = yp;
% 
%     for k = 1:max_iter
%         y_new = S(j) + (h/24)*...
%             (9*f(L(j+1),yc) ...
%             +19*f(L(j),S(j)) ...
%             -5*f(L(j-1),S(j-1)) ...
%             +f(L(j-2),S(j-2)));
% 
%         if abs(y_new - yc) < tol
%             break;
%         end
% 
%         yc = y_new;
%     end
% 
%     S(j+1) = yc;
% 
%     fprintf('j= %2d, Corrected Solution = %12.12f (iters=%d)\n\n', j, S(j+1), k);
% end
% 
% % %% ---- Exact solution ----
% % syms u(t)
% % usol(t) = dsolve(diff(u,t) == (t^2 + u(t)^2), u(0) == 1);
% % 
% % Error = abs(double(usol(b)) - S(N+1));
% % 
% % fprintf('Final Error = %12.12e\n\n', Error);
% 
% %% ---- solution using ode45 ----
% [x_ref, y_ref] = ode45(f, [a b], 1);
% %[x_ref, y_ref] = ode45(f, [a b], 1, ...
%     %odeset('RelTol',1e-12,'AbsTol',1e-14));
% y_exact = y_ref(end);
% 
% Error = abs(y_exact - S(N+1));
% 
% fprintf('Final Error = %12.12e\n\n', Error);




