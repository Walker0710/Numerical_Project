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

% Storage
x = zeros(N+1,1);  % Hares
y = zeros(N+1,1);  % Lynx

% Initial values (year 1900)
x(1) = 12.82;
y(1) = 7.13;

% Explicit Euler loop
for n = 1:N
    x(n+1) = x(n) + h*(alpha*x(n) - beta*x(n)*y(n));
    y(n+1) = y(n) + h*(delta*x(n)*y(n) - gamma*y(n));
end

% Sample Euler solution at each year (0,1,2,...,35)
steps_per_year = round(1/h);
year_idx = 1:steps_per_year:(N+1);

hares_euler_yearly = x(year_idx);
lynx_euler_yearly  = y(year_idx);

% Error metrics vs paper data
rmse_hares = sqrt(mean((hares_euler_yearly - hares_obs).^2));
rmse_lynx  = sqrt(mean((lynx_euler_yearly - lynx_obs).^2));

% Display results
fprintf('Explicit Euler (h = %.1f) using paper Eq. (12) parameters\n', h);
fprintf('RMSE Hares = %.4f\n', rmse_hares);
fprintf('RMSE Lynx  = %.4f\n\n', rmse_lynx);

fprintf(' Year   Hare_obs   Hare_euler   Lynx_obs   Lynx_euler\n');
fprintf('-------------------------------------------------------\n');
for i = 1:length(years)
    fprintf('%4d   %8.2f   %10.4f   %8.2f   %10.4f\n', ...
        years(i), hares_obs(i), hares_euler_yearly(i), lynx_obs(i), lynx_euler_yearly(i));
end
