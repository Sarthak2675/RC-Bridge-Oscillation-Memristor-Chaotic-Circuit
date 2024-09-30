% Parameters
a = 40;     % Adjust as needed
b = 0.15;   % Adjust as needed
c = -0.2;   % Adjust as needed
W_factor = 0.0614;
W_offset = 0.03;

% Time span and step reduction
tspan = linspace(0, 100, 500);  % Increased number of time points for better resolution

% Initial conditions for the two cases
initial_conditions_1 = [1; 0; 0; 0];  % First initial condition set
initial_conditions_2 = [-1; 0; 0; 0]; % Second initial condition set

% ODE function
ode_func = @(t, y, a) [
    a * (b * (3.2 - y(2) - y(3) - y(1)) - c * y(1) - y(1) * (W_offset + W_factor * abs(y(4))));
    (3.2 - 2) * y(2) - y(3);
    b * (3.2 * y(2) - y(1) - y(3)) + (3.2 - 1) * y(2) - y(3);
    y(1)
];

% Set ODE solver tolerances
opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'MaxStep', 1e-2);  % Solver options

%% Solve for both initial conditions
[t1, sol1] = ode15s(@(t, y) ode_func(t, y, a), tspan, initial_conditions_1, opts);
[t2, sol2] = ode15s(@(t, y) ode_func(t, y, a), tspan, initial_conditions_2, opts);

% Extract solutions for both initial conditions
x1 = sol1(:, 1); y1 = sol1(:, 2); z1 = sol1(:, 3); omega1 = sol1(:, 4);
x2 = sol2(:, 1); y2 = sol2(:, 2); z2 = sol2(:, 3); omega2 = sol2(:, 4);

%% Plotting Phase Portraits with different initial conditions
figure;
subplot(1, 2, 1);
plot(z1, x1, 'b'); hold on;
plot(z2, x2, 'r');
xlabel('z');
ylabel('x');
title('Phase Portrait (z vs x)');
legend('Initial cond 1', 'Initial cond 2');
grid on;

subplot(1, 2, 2);
plot(y1, x1, 'b'); hold on;
plot(y2, x2, 'r');
xlabel('y');
ylabel('x');
title('Phase Portrait (x vs y)');
legend('Initial cond 1', 'Initial cond 2');
grid on;

figure;
subplot(1, 2, 1);
plot(y1, z1, 'b'); hold on;
plot(y2, z2, 'r');
xlabel('y');
ylabel('z');
title('Phase Portrait (z vs y)');
legend('Initial cond 1', 'Initial cond 2');
grid on;

subplot(1, 2, 2);
plot(omega1, x1, 'b'); hold on;
plot(omega2, x2, 'r');
xlabel('omega');
ylabel('x');
title('Phase Portrait (omega vs x)');
legend('Initial cond 1', 'Initial cond 2');
grid on;

figure;
plot3(x1, y1, omega1, 'b'); hold on;
plot3(x2, y2, omega2, 'r');
xlabel('x');
ylabel('y');
zlabel('omega');
title('3D Phase Portrait (x, y, omega)');
legend('Initial cond 1', 'Initial cond 2');
grid on;

%% Poincaré Map for x = 0 with different initial conditions
crossings_x1 = find(diff(sign(x1)) ~= 0);
y_values_at_x0_1 = y1(crossings_x1);
omega_values_at_x0_1 = omega1(crossings_x1);

crossings_x2 = find(diff(sign(x2)) ~= 0);
y_values_at_x0_2 = y2(crossings_x2);
omega_values_at_x0_2 = omega2(crossings_x2);

figure;
plot(y_values_at_x0_1, omega_values_at_x0_1, 'bo'); hold on;
plot(y_values_at_x0_2, omega_values_at_x0_2, 'ro');
xlabel('y (at x crossing)');
ylabel('omega (at x crossing)');
title('Poincaré Map for x = 0');
legend('Initial cond 1', 'Initial cond 2');
grid on;

%% Poincaré Map for y = 0 with different initial conditions
crossings_y1 = find(diff(sign(y1)) ~= 0);
x_values_at_y0_1 = x1(crossings_y1);
omega_values_at_y0_1 = omega1(crossings_y1);

crossings_y2 = find(diff(sign(y2)) ~= 0);
x_values_at_y0_2 = x2(crossings_y2);
omega_values_at_y0_2 = omega2(crossings_y2);

figure;
plot(x_values_at_y0_1, omega_values_at_y0_1, 'bo'); hold on;
plot(x_values_at_y0_2, omega_values_at_y0_2, 'ro');
xlabel('x (at y crossing)');
ylabel('omega (at y crossing)');
title('Poincaré Map for y = 0');
legend('Initial cond 1', 'Initial cond 2');
grid on;

%% Bifurcation Diagram: Collect steady-state values for varying 'a' for both initial conditions
a_values = linspace(20, 60, 100);  % Vary 'a' from 20 to 60
bifurcation_values_1 = [];
bifurcation_values_2 = [];

for a = a_values
    % Solve for initial_conditions_1
    [~, sol1] = ode15s(@(t, y) ode_func(t, y, a), tspan, initial_conditions_1, opts);
    bifurcation_values_1 = [bifurcation_values_1; a * ones(length(tspan), 1), sol1(:, 1)];  % Store 'a' and x values for initial_conditions_1
    
    % Solve for initial_conditions_2
    [~, sol2] = ode15s(@(t, y) ode_func(t, y, a), tspan, initial_conditions_2, opts);
    bifurcation_values_2 = [bifurcation_values_2; a * ones(length(tspan), 1), sol2(:, 1)];  % Store 'a' and x values for initial_conditions_2
end

% Plot Bifurcation Diagram for both initial conditions
figure;
plot(bifurcation_values_1(:, 1), bifurcation_values_1(:, 2), 'b.'); hold on;
plot(bifurcation_values_2(:, 1), bifurcation_values_2(:, 2), 'r.');
xlabel('Parameter a');
ylabel('Steady-State x');
title('Bifurcation Diagram for varying a');
legend('Initial cond 1', 'Initial cond 2');
grid on;
