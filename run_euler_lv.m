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

% LV RHS (defined once and reused by both methods)
f1 = @(x,y) alpha*x - beta*x*y;   % dx/dt
f2 = @(x,y) delta*x*y - gamma*y;  % dy/dt

% Storage: Explicit Euler
x_euler = zeros(N+1,1);  % Hares
y_euler = zeros(N+1,1);  % Lynx

% Initial values (year 1900)
x_euler(1) = 12.82;
y_euler(1) = 7.13;

% Explicit Euler
for n = 1:N
    x_euler(n+1) = x_euler(n) + h*f1(x_euler(n), y_euler(n));
    y_euler(n+1) = y_euler(n) + h*f2(x_euler(n), y_euler(n));
end

% Storage: Modified Euler / Heun
x_heun = zeros(N+1,1);  % Hares
y_heun = zeros(N+1,1);  % Lynx

% Initial values (year 1900)
x_heun(1) = 12.82;
y_heun(1) = 7.13;

% Modified Euler / Heun
for n = 1:N
    x_tilde = x_heun(n) + h*f1(x_heun(n), y_heun(n));
    y_tilde = y_heun(n) + h*f2(x_heun(n), y_heun(n));

    x_heun(n+1) = x_heun(n) + (h/2)*(f1(x_heun(n), y_heun(n)) + f1(x_tilde, y_tilde));
    y_heun(n+1) = y_heun(n) + (h/2)*(f2(x_heun(n), y_heun(n)) + f2(x_tilde, y_tilde));
end

% Sample both methods at each year (0,1,2,...,35)
steps_per_year = round(1/h);
year_idx = 1:steps_per_year:(N+1);

hare_euler_yearly = x_euler(year_idx);
lynx_euler_yearly = y_euler(year_idx);
hare_heun_yearly  = x_heun(year_idx);
lynx_heun_yearly  = y_heun(year_idx);

% Display combined comparison table (paper -> Euler -> modified Euler)
fprintf('Lotka-Volterra Comparison (h = %.1f)\n', h);
fprintf('Order shown: paper value, Euler value, modified Euler value\n\n');

fprintf(' Year   Hare_paper   Hare_euler   Hare_modified_euler   Lynx_paper   Lynx_euler   Lynx_modified_euler\n');
fprintf('--------------------------------------------------------------------------------------------------------\n');
for i = 1:length(years)
    fprintf('%4d   %10.2f   %10.4f   %19.4f   %10.2f   %10.4f   %19.4f\n', ...
        years(i), hares_obs(i), hare_euler_yearly(i), hare_heun_yearly(i), ...
        lynx_obs(i), lynx_euler_yearly(i), lynx_heun_yearly(i));
end
