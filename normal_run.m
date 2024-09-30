% Parameters
a = 40;     % Adjust as needed
b = 0.15;   % Adjust as needed
c = -0.2;   % Adjust as needed
W_factor = 0.0614;
W_offset = 0.03;

% Time span and step reduction
tspan = linspace(0, 100, 500);  % Increased number of time points for better resolution

% Initial conditions
initial_conditions = [0.01; 0; 0; 0];  % Initial conditions for [x; y; z; omega]

% ODE function
ode_func = @(t, y, A) [
    a * (b * (A - y(2) - y(3) - y(1)) - c * y(1) - y(1) * (W_offset + W_factor * abs(y(4))));
    (A - 2) * y(2) - y(3);
    b * (A * y(2) - y(1) - y(3)) + (A - 1) * y(2) - y(3);
    y(1)
];

% Set ODE solver tolerances
opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'MaxStep', 1e-2);  % Further relaxed tolerances and increased max step size

% Solve the system of ODEs using ode15s for stiff problems
tic;
[t, sol] = ode15s(@(t, y) ode_func(t, y, 3.2), tspan, initial_conditions, opts);  % Fixed A for initial solve
elapsed_time_solving = toc;
disp(['Time to solve ODEs: ' num2str(elapsed_time_solving) ' seconds']);

% Extract solutions
x = sol(:, 1);
y = sol(:, 2);
z = sol(:, 3);
omega = sol(:, 4);

%% Plotting Phase Portraits
figure;
subplot(1, 2, 1);
plot(z, x);
xlabel('x');
ylabel('z');
title('Phase Portrait (z vs x)');
grid on;

subplot(1, 2, 2);
plot(y, x);
xlabel('y');
ylabel('x');
title('Phase Portrait (x vs y)');
grid on;

figure;
subplot(1, 2, 1);
plot(y, z);
xlabel('y');
ylabel('z');
title('Phase Portrait (z vs y)');
grid on;

subplot(1, 2, 2);
plot(omega, x);
xlabel('omega');
ylabel('x');
title('Phase Portrait (omega vs x)');
grid on;

figure;
plot3(x, y, omega);
xlabel('x');
ylabel('y');
zlabel('omega');
title('3D Phase Portrait (x, y, omega)');
grid on;

%% Poincaré Map
% Find x crossings
crossings_x = find(diff(sign(x)) ~= 0);
y_values_at_x0 = y(crossings_x);
omega_values_at_x0 = omega(crossings_x);

figure;
plot(y_values_at_x0, omega_values_at_x0, 'o');
xlabel('y (at x crossing)');
ylabel('omega (at x crossing)');
title('Poincaré Map for x = 0');
grid on;

%% Poincaré Map for y = 0
% Find y crossings
crossings_y = find(diff(sign(y)) ~= 0);
x_values_at_y0 = x(crossings_y);
omega_values_at_y0 = omega(crossings_y);

% Plot Poincaré Map for y = 0
figure;
plot(x_values_at_y0, omega_values_at_y0, 'o');
xlabel('x (at y crossing)');
ylabel('omega (at y crossing)');
title('Poincaré Map for y = 0');
grid on;

%% Lyapunov Exponent Calculation Optimization
tic;  % Start timer for Lyapunov exponent calculation
N = length(t);
lyapunov_exponents = zeros(1, 3);
step_size = 10;  % Calculate every 10 steps to optimize memory and computation

for i = 1:step_size:N-step_size
    % Calculate divergence in each direction
    divergence_x = log(abs(diff(x(i:i+step_size))));
    divergence_y = log(abs(diff(y(i:i+step_size))));
    divergence_z = log(abs(diff(z(i:i+step_size))));
    
    if ~isnan(divergence_x)
        lyapunov_exponents(1) = lyapunov_exponents(1) + sum(divergence_x);
    end
    if ~isnan(divergence_y)
        lyapunov_exponents(2) = lyapunov_exponents(2) + sum(divergence_y);
    end
    if ~isnan(divergence_z)
        lyapunov_exponents(3) = lyapunov_exponents(3) + sum(divergence_z);
    end
end

% Average Lyapunov Exponents
lyapunov_exponents = lyapunov_exponents / (N/step_size);
elapsed_time_lyapunov = toc;
disp(['Time to calculate Lyapunov Exponents: ' num2str(elapsed_time_lyapunov) ' seconds']);

disp('Lyapunov Exponents:');
disp(lyapunov_exponents);

% Calculate Lyapunov Spectrum over a range of A values
A_values = linspace(2, 4, 100);  % Vary A from 2 to 4
lyapunov_spectrum = zeros(length(A_values), 3);

for i = 1:length(A_values)
    [t, sol] = ode15s(@(t, y) ode_func(t, y, A_values(i)), tspan, initial_conditions, opts);
    x_sol = sol(:, 1);
    y_sol = sol(:, 2);
    z_sol = sol(:, 3);
    
    % Calculate Lyapunov Exponents for each A
    divergence_x = log(abs(diff(x_sol)));
    divergence_y = log(abs(diff(y_sol)));
    divergence_z = log(abs(diff(z_sol)));
    
    lyapunov_spectrum(i, 1) = sum(divergence_x) / (length(x_sol) - 1);
    lyapunov_spectrum(i, 2) = sum(divergence_y) / (length(y_sol) - 1);
    lyapunov_spectrum(i, 3) = sum(divergence_z) / (length(z_sol) - 1);
end

% Plot Lyapunov Spectrum
figure;
plot(A_values, lyapunov_spectrum);
xlabel('Parameter A');
ylabel('Lyapunov Exponents');
title('Lyapunov Spectrum');
legend('L_1', 'L_2', 'L_3');
grid on;

% Bifurcation Diagram: Collect steady-state values for varying A
bifurcation_values = [];

for A = A_values
    [t, sol] = ode15s(@(t, y) ode_func(t, y, A), tspan, initial_conditions, opts);
    bifurcation_values = [bifurcation_values; A * ones(length(t), 1), sol(:, 1)];  % Store A and x values
end

% Plot Bifurcation Diagram
figure;
plot(bifurcation_values(:, 1), bifurcation_values(:, 2), '.');
xlabel('Parameter A');
ylabel('Steady-State x');
title('Bifurcation Diagram');
grid on;
