% Parameters
b = 0.15;   % Adjust as needed
c = -0.2;   % Adjust as needed
W_factor = 0.0614;
W_offset = 0.03;

% Time span and step reduction
tspan = linspace(0, 100, 500);  % Increased number of time points for better resolution

% Initial conditions
initial_conditions = [0.01; 0; 0; 0];  % Initial conditions for [x; y; z; omega]

% ODE function with varying 'a'
ode_func = @(t, y, a) [
    a * (b * (3.2 - y(2) - y(3) - y(1)) - c * y(1) - y(1) * (W_offset + W_factor * abs(y(4))));  % A fixed, a varied
    (3.2 - 2) * y(2) - y(3);
    b * (3.2 * y(2) - y(1) - y(3)) + (3.2 - 1) * y(2) - y(3);
    y(1)
];

% Set ODE solver tolerances
opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'MaxStep', 1e-2);  % Further relaxed tolerances and increased max step size

%% Calculate Lyapunov Spectrum over a range of 'a' values
a_values = linspace(20, 60, 100);  % Vary 'a' from 20 to 60
lyapunov_spectrum = zeros(length(a_values), 3);

for i = 1:length(a_values)
    [t, sol] = ode15s(@(t, y) ode_func(t, y, a_values(i)), tspan, initial_conditions, opts);
    x_sol = sol(:, 1);
    y_sol = sol(:, 2);
    z_sol = sol(:, 3);
    
    % Calculate Lyapunov Exponents for each 'a'
    divergence_x = log(abs(diff(x_sol)));
    divergence_y = log(abs(diff(y_sol)));
    divergence_z = log(abs(diff(z_sol)));
    
    lyapunov_spectrum(i, 1) = sum(divergence_x) / (length(x_sol) - 1);
    lyapunov_spectrum(i, 2) = sum(divergence_y) / (length(y_sol) - 1);
    lyapunov_spectrum(i, 3) = sum(divergence_z) / (length(z_sol) - 1);
end

% Plot Lyapunov Spectrum for varying 'a'
figure;
plot(a_values, lyapunov_spectrum);
xlabel('Parameter a');
ylabel('Lyapunov Exponents');
title('Lyapunov Spectrum for varying a');
legend('L_1', 'L_2', 'L_3');
grid on;

%% Bifurcation Diagram: Collect steady-state values for varying 'a'
bifurcation_values = [];

for a = a_values
    [t, sol] = ode15s(@(t, y) ode_func(t, y, a), tspan, initial_conditions, opts);
    bifurcation_values = [bifurcation_values; a * ones(length(t), 1), sol(:, 1)];  % Store 'a' and x values
end

% Plot Bifurcation Diagram for varying 'a'
figure;
plot(bifurcation_values(:, 1), bifurcation_values(:, 2), '.');
xlabel('Parameter a');
ylabel('Steady-State x');
title('Bifurcation Diagram for varying a');
grid on;

%% Solve ODE system for a specific 'a' (fixed) to plot phase portraits, Poincaré maps, etc.
fixed_a = 40;  % You can choose a specific value for 'a'
[t, sol] = ode15s(@(t, y) ode_func(t, y, fixed_a), tspan, initial_conditions, opts);

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
title(['Phase Portrait (z vs x) for a = ' num2str(fixed_a)]);
grid on;

subplot(1, 2, 2);
plot(y, x);
xlabel('y');
ylabel('x');
title(['Phase Portrait (x vs y) for a = ' num2str(fixed_a)]);
grid on;

figure;
subplot(1, 2, 1);
plot(y, z);
xlabel('y');
ylabel('z');
title(['Phase Portrait (z vs y) for a = ' num2str(fixed_a)]);
grid on;

subplot(1, 2, 2);
plot(omega, x);
xlabel('omega');
ylabel('x');
title(['Phase Portrait (omega vs x) for a = ' num2str(fixed_a)]);
grid on;

figure;
plot3(x, y, omega);
xlabel('x');
ylabel('y');
zlabel('omega');
title(['3D Phase Portrait (x, y, omega) for a = ' num2str(fixed_a)]);
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
title(['Poincaré Map for x = 0 and a = ' num2str(fixed_a)]);
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
title(['Poincaré Map for y = 0 and a = ' num2str(fixed_a)]);
grid on;
