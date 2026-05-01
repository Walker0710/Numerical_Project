clear; clc;

% Explicit Euler for the Lotka-Volterra system
% dx/dt = alpha*x - beta*x*y
% dy/dt = delta*x*y - gamma*y

% Parameters from paper Eq. (12)
alpha = 0.6923;
beta  = 0.0333;
gamma = 0.7627;
delta = 0.02684;

% Observed paper data (1900 to 1935)
years = (1900:1935)';
hares_obs = [ ...
    12.82; 4.72; 4.73; 37.22; 69.72; 57.78; 28.68; 23.37; 21.54; 26.34; ...
    53.10; 68.48; 75.58; 57.92; 40.97; 24.95; 12.59; 4.97; 4.50; 11.21; ...
    56.60; 69.63; 77.74; 80.53; 73.38; 36.93; 4.64; 2.54; 1.80; 2.39; ...
    4.23; 19.52; 82.11; 89.76; 81.66; 15.76 ...
];
lynx_obs = [ ...
    7.13; 9.47; 14.86; 31.47; 60.57; 63.51; 54.70; 6.30; 3.41; 5.44; ...
    11.65; 20.35; 32.88; 39.55; 43.36; 40.83; 30.36; 17.18; 6.82; 3.19; ...
    3.52; 9.94; 20.30; 31.99; 42.36; 49.08; 53.99; 52.25; 37.70; 19.14; ...
    6.98; 8.31; 16.01; 24.82; 29.70; 35.40 ...
];

% Time setup: small step for stability, then sample at each year
h = 0.1;
total_years = years(end) - years(1);
N = round(total_years/h);

% LV RHS (defined once and reused by all methods)
f1 = @(t,y1,y2) alpha*y1 - beta*y1*y2;   % dx/dt
f2 = @(t,y1,y2) delta*y1*y2 - gamma*y2;  % dy/dt

% Storage: Explicit Euler
x_euler = zeros(N+1,1);  % Hares
y_euler = zeros(N+1,1);  % Lynx

% Initial values (year 1900)
x_euler(1) = 12.82;
y_euler(1) = 7.13;

% Explicit Euler
for n = 1:N
    t_n = (n-1)*h;
    x_euler(n+1) = x_euler(n) + h*f1(t_n, x_euler(n), y_euler(n));
    y_euler(n+1) = y_euler(n) + h*f2(t_n, x_euler(n), y_euler(n));
end

% Storage: Modified Euler / Heun
x_heun = zeros(N+1,1);  % Hares
y_heun = zeros(N+1,1);  % Lynx

% Initial values (year 1900)
x_heun(1) = 12.82;
y_heun(1) = 7.13;

% Modified Euler / Heun
for n = 1:N
    t_n = (n-1)*h;
    x_tilde = x_heun(n) + h*f1(t_n, x_heun(n), y_heun(n));
    y_tilde = y_heun(n) + h*f2(t_n, x_heun(n), y_heun(n));

    x_heun(n+1) = x_heun(n) + (h/2)*(f1(t_n, x_heun(n), y_heun(n)) + f1(t_n+h, x_tilde, y_tilde));
    y_heun(n+1) = y_heun(n) + (h/2)*(f2(t_n, x_heun(n), y_heun(n)) + f2(t_n+h, x_tilde, y_tilde));
end

% Storage: RK4
x_rk4 = zeros(N+1,1);  % Hares
y_rk4 = zeros(N+1,1);  % Lynx

% Initial values (year 1900)
x_rk4(1) = 12.82;
y_rk4(1) = 7.13;

% RK4
for n = 1:N
    t_n = (n-1)*h;
    y1 = x_rk4(n);
    y2 = y_rk4(n);

    k1_1 = h * f1(t_n, y1, y2);
    k1_2 = h * f2(t_n, y1, y2);

    k2_1 = h * f1(t_n + h/2, y1 + k1_1/2, y2 + k1_2/2);
    k2_2 = h * f2(t_n + h/2, y1 + k1_1/2, y2 + k1_2/2);

    k3_1 = h * f1(t_n + h/2, y1 + k2_1/2, y2 + k2_2/2);
    k3_2 = h * f2(t_n + h/2, y1 + k2_1/2, y2 + k2_2/2);

    k4_1 = h * f1(t_n + h, y1 + k3_1, y2 + k3_2);
    k4_2 = h * f2(t_n + h, y1 + k3_1, y2 + k3_2);

    x_rk4(n+1) = y1 + (k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6;
    y_rk4(n+1) = y2 + (k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6;
end

% Storage: AB4-AM4 Predictor-Corrector (single corrector)
x_abm = zeros(N+1,1);  % Hares
y_abm = zeros(N+1,1);  % Lynx

% Initial values (year 1900)
x_abm(1) = 12.82;
y_abm(1) = 7.13;

% Startup: RK4 for first 3 steps
for n = 1:3
    t_n = (n-1)*h;
    y1 = x_abm(n);
    y2 = y_abm(n);

    k1_1 = h * f1(t_n, y1, y2);
    k1_2 = h * f2(t_n, y1, y2);

    k2_1 = h * f1(t_n + h/2, y1 + k1_1/2, y2 + k1_2/2);
    k2_2 = h * f2(t_n + h/2, y1 + k1_1/2, y2 + k1_2/2);

    k3_1 = h * f1(t_n + h/2, y1 + k2_1/2, y2 + k2_2/2);
    k3_2 = h * f2(t_n + h/2, y1 + k2_1/2, y2 + k2_2/2);

    k4_1 = h * f1(t_n + h, y1 + k3_1, y2 + k3_2);
    k4_2 = h * f2(t_n + h, y1 + k3_1, y2 + k3_2);

    x_abm(n+1) = y1 + (k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6;
    y_abm(n+1) = y2 + (k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6;
end

% AB4 predictor + AM4 corrector (single pass)
for n = 4:N
    t_n = (n-1)*h;
    t_nm1 = (n-2)*h;
    t_nm2 = (n-3)*h;
    t_nm3 = (n-4)*h;
    t_np1 = n*h;

    f1_n   = f1(t_n,   x_abm(n),   y_abm(n));
    f1_nm1 = f1(t_nm1, x_abm(n-1), y_abm(n-1));
    f1_nm2 = f1(t_nm2, x_abm(n-2), y_abm(n-2));
    f1_nm3 = f1(t_nm3, x_abm(n-3), y_abm(n-3));

    f2_n   = f2(t_n,   x_abm(n),   y_abm(n));
    f2_nm1 = f2(t_nm1, x_abm(n-1), y_abm(n-1));
    f2_nm2 = f2(t_nm2, x_abm(n-2), y_abm(n-2));
    f2_nm3 = f2(t_nm3, x_abm(n-3), y_abm(n-3));

    % Predictor (AB4)
    x_pred = x_abm(n) + (h/24)*(55*f1_n - 59*f1_nm1 + 37*f1_nm2 - 9*f1_nm3);
    y_pred = y_abm(n) + (h/24)*(55*f2_n - 59*f2_nm1 + 37*f2_nm2 - 9*f2_nm3);

    % Corrector (AM4, single pass)
    x_abm(n+1) = x_abm(n) + (h/24)*(9*f1(t_np1, x_pred, y_pred) + 19*f1_n - 5*f1_nm1 + f1_nm2);
    y_abm(n+1) = y_abm(n) + (h/24)*(9*f2(t_np1, x_pred, y_pred) + 19*f2_n - 5*f2_nm1 + f2_nm2);
end

% Sample both methods at each year (0,1,2,...,35)
steps_per_year = round(1/h);
year_idx = 1:steps_per_year:(N+1);

hare_euler_yearly = x_euler(year_idx);
lynx_euler_yearly = y_euler(year_idx);
hare_heun_yearly  = x_heun(year_idx);
lynx_heun_yearly  = y_heun(year_idx);
hare_rk4_yearly   = x_rk4(year_idx);
lynx_rk4_yearly   = y_rk4(year_idx);
hare_abm_yearly   = x_abm(year_idx);
lynx_abm_yearly   = y_abm(year_idx);

% Storage: ODE45
f_ode45 = @(t, Y) [alpha*Y(1) - beta*Y(1)*Y(2); delta*Y(1)*Y(2) - gamma*Y(2)];
[~, Y_ode45] = ode45(f_ode45, years, [12.82; 7.13]);
hare_ode45_yearly = Y_ode45(:,1);
lynx_ode45_yearly = Y_ode45(:,2);

% RMSE for each method
rmse_hare_euler = sqrt(mean((hare_euler_yearly - hares_obs).^2));
rmse_lynx_euler = sqrt(mean((lynx_euler_yearly - lynx_obs).^2));

rmse_hare_heun = sqrt(mean((hare_heun_yearly - hares_obs).^2));
rmse_lynx_heun = sqrt(mean((lynx_heun_yearly - lynx_obs).^2));

rmse_hare_rk4 = sqrt(mean((hare_rk4_yearly - hares_obs).^2));
rmse_lynx_rk4 = sqrt(mean((lynx_rk4_yearly - lynx_obs).^2));

rmse_hare_abm = sqrt(mean((hare_abm_yearly - hares_obs).^2));
rmse_lynx_abm = sqrt(mean((lynx_abm_yearly - lynx_obs).^2));

rmse_hare_ode45 = sqrt(mean((hare_ode45_yearly - hares_obs).^2));
rmse_lynx_ode45 = sqrt(mean((lynx_ode45_yearly - lynx_obs).^2));

% Display combined comparison table
fprintf('Lotka-Volterra Comparison (h = %.1f)\n', h);
fprintf('Order shown: paper value, Euler value, modified Euler value, RK4 value, ABM value, ODE45 value\n\n');

fprintf(' Year   Hare_paper   Hare_euler   Hare_modified_euler   Hare_rk4   Hare_abm   Hare_ode45   Lynx_paper   Lynx_euler   Lynx_modified_euler   Lynx_rk4   Lynx_abm   Lynx_ode45\n');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
for i = 1:length(years)
    fprintf('%4d   %10.2f   %10.4f   %19.4f   %9.4f   %8.4f   %10.4f   %10.2f   %10.4f   %19.4f   %9.4f   %8.4f   %10.4f\n', ...
        years(i), hares_obs(i), hare_euler_yearly(i), hare_heun_yearly(i), hare_rk4_yearly(i), hare_abm_yearly(i), hare_ode45_yearly(i), ...
        lynx_obs(i), lynx_euler_yearly(i), lynx_heun_yearly(i), lynx_rk4_yearly(i), lynx_abm_yearly(i), lynx_ode45_yearly(i));
end

% Separate RMSE comparison table
fprintf('\nRMSE Comparison Table\n');
fprintf('---------------------------------------------\n');
fprintf(' Method                 RMSE_Hare   RMSE_Lynx\n');
fprintf('---------------------------------------------\n');
fprintf(' Euler                  %9.4f   %9.4f\n', rmse_hare_euler, rmse_lynx_euler);
fprintf(' Modified Euler (Heun)  %9.4f   %9.4f\n', rmse_hare_heun, rmse_lynx_heun);
fprintf(' RK4                    %9.4f   %9.4f\n', rmse_hare_rk4, rmse_lynx_rk4);
fprintf(' AB4-AM4                %9.4f   %9.4f\n', rmse_hare_abm, rmse_lynx_abm);
fprintf(' ODE45                  %9.4f   %9.4f\n', rmse_hare_ode45, rmse_lynx_ode45);

%% Plotting

% 1. Paper vs Euler
figure('Name', 'Paper vs Euler', 'NumberTitle', 'off');
subplot(2,1,1);
plot(years, hares_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, hare_euler_yearly, 'b-*', 'LineWidth', 1.5);
title('Hare Population: Paper vs Euler');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed (Paper)', 'Euler');
grid on;

subplot(2,1,2);
plot(years, lynx_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, lynx_euler_yearly, 'r-*', 'LineWidth', 1.5);
title('Lynx Population: Paper vs Euler');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed (Paper)', 'Euler');
grid on;
saveas(gcf, 'plot1_paper_vs_euler.png');

% 2. Paper vs Modified Euler (Heun)
figure('Name', 'Paper vs Modified Euler', 'NumberTitle', 'off');
subplot(2,1,1);
plot(years, hares_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, hare_heun_yearly, 'b-s', 'LineWidth', 1.5);
title('Hare Population: Paper vs Modified Euler');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed (Paper)', 'Modified Euler');
grid on;

subplot(2,1,2);
plot(years, lynx_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, lynx_heun_yearly, 'r-s', 'LineWidth', 1.5);
title('Lynx Population: Paper vs Modified Euler');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed (Paper)', 'Modified Euler');
grid on;
saveas(gcf, 'plot2_paper_vs_modified_euler.png');

% 3. Paper vs RK4
figure('Name', 'Paper vs RK4', 'NumberTitle', 'off');
subplot(2,1,1);
plot(years, hares_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, hare_rk4_yearly, 'b-^', 'LineWidth', 1.5);
title('Hare Population: Paper vs RK4');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed (Paper)', 'RK4');
grid on;

subplot(2,1,2);
plot(years, lynx_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, lynx_rk4_yearly, 'r-^', 'LineWidth', 1.5);
title('Lynx Population: Paper vs RK4');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed (Paper)', 'RK4');
grid on;
saveas(gcf, 'plot3_paper_vs_rk4.png');

% 4. Paper vs AB4-AM4
figure('Name', 'Paper vs AB4-AM4', 'NumberTitle', 'off');
subplot(2,1,1);
plot(years, hares_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, hare_abm_yearly, 'b-d', 'LineWidth', 1.5);
title('Hare Population: Paper vs AB4-AM4');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed (Paper)', 'AB4-AM4');
grid on;

subplot(2,1,2);
plot(years, lynx_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, lynx_abm_yearly, 'r-d', 'LineWidth', 1.5);
title('Lynx Population: Paper vs AB4-AM4');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed (Paper)', 'AB4-AM4');
grid on;
saveas(gcf, 'plot4_paper_vs_ab4_am4.png');

% 5. Paper vs ODE45
figure('Name', 'Paper vs ODE45', 'NumberTitle', 'off');
subplot(2,1,1);
plot(years, hares_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, hare_ode45_yearly, 'b-p', 'LineWidth', 1.5);
title('Hare Population: Paper vs ODE45');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed (Paper)', 'ODE45');
grid on;

subplot(2,1,2);
plot(years, lynx_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, lynx_ode45_yearly, 'r-p', 'LineWidth', 1.5);
title('Lynx Population: Paper vs ODE45');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed (Paper)', 'ODE45');
grid on;
saveas(gcf, 'plot5_paper_vs_ode45.png');

% 6. Paper vs All Methods
figure('Name', 'Paper vs All Methods', 'NumberTitle', 'off');
subplot(2,1,1);
plot(years, hares_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, hare_euler_yearly, 'b-*', 'LineWidth', 1.5);
plot(years, hare_heun_yearly, 'g-s', 'LineWidth', 1.5);
plot(years, hare_rk4_yearly, 'm-^', 'LineWidth', 1.5);
plot(years, hare_abm_yearly, 'c-d', 'LineWidth', 1.5);
plot(years, hare_ode45_yearly, 'y-p', 'LineWidth', 1.5);
title('Hare Population: All Methods vs Paper');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed', 'Euler', 'Modified Euler', 'RK4', 'AB4-AM4', 'ODE45', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(years, lynx_obs, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
plot(years, lynx_euler_yearly, 'b-*', 'LineWidth', 1.5);
plot(years, lynx_heun_yearly, 'g-s', 'LineWidth', 1.5);
plot(years, lynx_rk4_yearly, 'm-^', 'LineWidth', 1.5);
plot(years, lynx_abm_yearly, 'c-d', 'LineWidth', 1.5);
plot(years, lynx_ode45_yearly, 'y-p', 'LineWidth', 1.5);
title('Lynx Population: All Methods vs Paper');
xlabel('Year'); ylabel('Population (thousands)');
legend('Observed', 'Euler', 'Modified Euler', 'RK4', 'AB4-AM4', 'ODE45', 'Location', 'best');
grid on;
saveas(gcf, 'plot6_paper_vs_all_methods.png');
