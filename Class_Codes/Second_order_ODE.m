clear; clc; close all;

% PROBLEM:
% y'' + 4y = cos(t),  y(0)=1,  y'(0)=0
% Let y1=y and y2=y'; Then y1'=y2 and y2'=cost - 4*y1;

% PARAMETERS
t0 = 0;
tf = 0.4;
h  = 0.2;

N = round((tf - t0)/h);

% INITIAL CONDITIONS
y1 = 1;    % y(0)
y2 = 0;    % y'(0)
t  = t0;

% STORAGE (RK4)
T_RK  = zeros(N+1,1);
Y1_RK = zeros(N+1,1);
Y2_RK = zeros(N+1,1);

T_RK(1)  = t;
Y1_RK(1) = y1;
Y2_RK(1) = y2;

% DEFINE SYSTEM
f1 = @(t,y1,y2) y2;
f2 = @(t,y1,y2) cos(t) - 4*y1;

% RK4 METHOD
for i = 1:N
    
    k1_1 = h * f1(t, y1, y2);
    k1_2 = h * f2(t, y1, y2);
    
    k2_1 = h * f1(t + h/2, y1 + k1_1/2, y2 + k1_2/2);
    k2_2 = h * f2(t + h/2, y1 + k1_1/2, y2 + k1_2/2);
    
    k3_1 = h * f1(t + h/2, y1 + k2_1/2, y2 + k2_2/2);
    k3_2 = h * f2(t + h/2, y1 + k2_1/2, y2 + k2_2/2);
    
    k4_1 = h * f1(t + h, y1 + k3_1, y2 + k3_2);
    k4_2 = h * f2(t + h, y1 + k3_1, y2 + k3_2);
    
    y1 = y1 + (k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6;
    y2 = y2 + (k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6;
    
    t = t + h;
    
    T_RK(i+1)  = t;
    Y1_RK(i+1) = y1;
    Y2_RK(i+1) = y2;
end

% ODE45 SOLUTION
% tspan = [t0 tf];
% y0 = [1; 0];
% 
% [t_ode, y_ode] = ode45(@ode_system, tspan, y0);

[t_ode,y_ode] = ode45(@(t,y) [y(2); cos(t)-4*y(1)], [0 0.4], [1;0]);

% DISPLAY RESULTS
disp('---------------- RK4 RESULT ----------------')
disp('   t        y(t)        y''(t)')
disp([T_RK Y1_RK Y2_RK])

disp('---------------- ODE45 RESULT ----------------')
disp('   t        y(t)        y''(t)')
disp([t_ode y_ode])

% FINAL VALUE COMPARISON
fprintf('\nAt t = %.2f:\n', tf);
fprintf('RK4   y = %.6f\n', Y1_RK(end));
fprintf('ODE45 y = %.6f\n', y_ode(end,1));

% PLOT COMPARISON
figure;
plot(T_RK, Y1_RK, 'ro-','LineWidth',1.5); hold on;
plot(t_ode, y_ode(:,1), 'b*-','LineWidth',1.5);

xlabel('t');
ylabel('y(t)');
title('Comparison: RK4 vs ODE45');
legend('RK4','ODE45');
grid on;

% FUNCTION FOR ODE45
function dydt = ode_system(t, y)

    dydt = zeros(2,1);

    % y(1) = y
    % y(2) = y'
    
    dydt(1) = y(2);                % y' = y2
    dydt(2) = cos(t) - 4*y(1);     % y'' = cos(t) - 4y

end