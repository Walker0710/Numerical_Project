clc; clear; close all;

% PROBLEM
% y'' = 6y^2 - x

f = @(x,y,z) 6*y.^2 - x;

a = 0; b = 1;
alpha = 1; beta = 5;

h =1/3;
x = a:h:b;

% SECANT METHOD (TWO SLOPES)

s1 = 1.0;
s2 = 2.0;

for iter = 1:10

    y1 = euler_solver(s1, x, h, alpha, f);
    F1 = y1(end) - beta;

    y2 = euler_solver(s2, x, h, alpha, f);
    F2 = y2(end) - beta;

    % Secant update
    s_new = s2 - F2*(s2 - s1)/(F2 - F1);

    if abs(s_new - s2) < 1e-8
        break;
    end

    s1 = s2;
    s2 = s_new;
end


% FINAL SOLUTION

y = euler_solver(s2, x, h, alpha, f);

fprintf('Final slope s = %.6f\n', s2);
fprintf('y(1/3) = %.6f\n', y(2));
fprintf('y(2/3) = %.6f\n', y(3));

disp('   x        y(x)')
disp([x' y'])


% EULER METHOD FUNCTION

function y = euler_solver(s, x, h, alpha, f)

N = length(x);

y = zeros(1,N);
z = zeros(1,N);

% Initial conditions
y(1) = alpha;
z(1) = s;

for i = 1:N-1
    xi = x(i);

    % Euler method
    y(i+1) = y(i) + h*z(i);
    z(i+1) = z(i) + h*f(xi, y(i), z(i));
end
end



% %%%% When we consider small step h<1/3 %%%%%%%%%
% 
% clc; clear; close all;
% 
% %% =========================================
% % PROBLEM
% % y'' = 6y^2 - x
% %% =========================================
% 
% f = @(x,y,z) 6*y.^2 - x;
% 
% a = 0; b = 1;
% alpha = 1; beta = 5;
% 
% h =0.1;
% x = a:h:b;
% 
% %% =========================================
% % SECANT METHOD (TWO SLOPES)
% %% =========================================
% 
% s1 = 1.0;
% s2 = 2.0;
% 
% for iter = 1:10
% 
%     y1 = euler_solver(s1, x, h, alpha, f);
%     F1 = y1(end) - beta;
% 
%     y2 = euler_solver(s2, x, h, alpha, f);
%     F2 = y2(end) - beta;
% 
%     % Secant update
%     s_new = s2 - F2*(s2 - s1)/(F2 - F1);
% 
%     if abs(s_new - s2) < 1e-8
%         break;
%     end
% 
%     s1 = s2;
%     s2 = s_new;
% end
% 
% % %% =========================================
% % % FINAL SOLUTION
% % %% =========================================
% % Correct evaluation using interpolation
% y = euler_solver(s2, x, h, alpha, f);
% y_13 = interp1(x, y, 1/3);
% y_23 = interp1(x, y, 2/3);
% 
% fprintf('Final slope s = %.6f\n', s2);
% fprintf('y(1/3) = %.6f\n', y_13);
% fprintf('y(2/3) = %.6f\n', y_23);
% 
% disp('   x        y(x)')
% disp([x' y'])
% 
% 
% %% =========================================
% % EULER METHOD FUNCTION
% %% =========================================
% 
% function y = euler_solver(s, x, h, alpha, f)
% 
% N = length(x);
% 
% y = zeros(1,N);
% z = zeros(1,N);
% 
% % Initial conditions
% y(1) = alpha;
% z(1) = s;
% 
% for i = 1:N-1
%     xi = x(i);
% 
%     % Euler method
%     y(i+1) = y(i) + h*z(i);
%     z(i+1) = z(i) + h*f(xi, y(i), z(i));
% end
% end