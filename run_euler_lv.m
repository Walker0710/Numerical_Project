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

% Display combined comparison table
fprintf('Lotka-Volterra Comparison (h = %.1f)\n', h);
fprintf('Order shown: paper value, Euler value, modified Euler value, RK4 value, ABM value\n\n');

fprintf(' Year   Hare_paper   Hare_euler   Hare_modified_euler   Hare_rk4   Hare_abm   Lynx_paper   Lynx_euler   Lynx_modified_euler   Lynx_rk4   Lynx_abm\n');
fprintf('----------------------------------------------------------------------------------------------------------------------------------------------------------\n');
for i = 1:length(years)
    fprintf('%4d   %10.2f   %10.4f   %19.4f   %9.4f   %8.4f   %10.2f   %10.4f   %19.4f   %9.4f   %8.4f\n', ...
        years(i), hares_obs(i), hare_euler_yearly(i), hare_heun_yearly(i), hare_rk4_yearly(i), hare_abm_yearly(i), ...
        lynx_obs(i), lynx_euler_yearly(i), lynx_heun_yearly(i), lynx_rk4_yearly(i), lynx_abm_yearly(i));
end
