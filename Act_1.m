% Parametros
Ip = 0.0079;  % Pendulum Moment of Inertia (kg*m^2)
Mc = 0.7031;  % Lumped Mass of the cart system (kg)
Lp = 0.3312;  % Distance from Pivot to Center of Gravity (m)
Mp = 0.32;    % Pendulum Mass with Fitting (kg)
Bc = 1.3;     % Equivalent Viscous Damping Coefficient at Motor Pinion (Ns/m)
q = 9.81;     % Gravitational Constant (m/s^2)
Bp = 0.0121;  % Viscous Damping Coefficient at Pendulum Axis (Nms/rad)

% Condiciones iniciales
alpha0 = deg2rad(1);  % Ángulo Inicial (rad)
x0 = 0;               % Posición inicial (m)
xdot0 = 0;            % Velocidad inicial (m/s)
alphadot0 = 0;        % Velocidad angular del pendulo (rad/s)

initial_conditions = [x0; xdot0; alpha0; alphadot0];

% Tiempo de simulación
tspan = [0, 10];  

inverted_pendulum = @(t, y) pendulum_dynamics(t, y, Ip, Mc, Mp, Lp, Bc, q, Bp);

[t, y] = ode45(inverted_pendulum, tspan, initial_conditions);

x = y(:, 1);       % Posición (m)
alpha = y(:, 3);   % Ángulo (rad)

figure;
set(gcf, 'Color', 'w'); 

plot(t, x, 'c-', 'LineWidth', 2); 
hold on;
plot(t, alpha, 'm-', 'LineWidth', 2); 
hold off;

xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Radianes', 'FontSize', 12, 'FontWeight', 'bold');
title('Respuesta combinada', 'FontSize', 14, 'FontWeight', 'bold');

legend('Posición del carro (x)', 'Ángulo del péndulo (\alpha)', 'Location', 'best');

grid on;

set(gca, 'FontSize', 12, 'FontWeight', 'bold');

figure;

subplot(2, 1, 1);
plot(t, x, 'b-', 'LineWidth', 2);  
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Posición del carro (m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Respuesta de la posición del carro', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% Plot pendulum angle
subplot(2, 1, 2);
plot(t, alpha, 'r-', 'LineWidth', 2); 
xlabel('Tiempo (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Ángulo del péndulo (rad)', 'FontSize', 12, 'FontWeight', 'bold');
title('Respuesta del ángulo del péndulo', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% Definicion de las funciones
function dydt = pendulum_dynamics(~, y, Ip, Mc, Mp, Lp, Bc, q, Bp)
    % Extraer variables de estado
    x = y(1);          % Posición del carro (m)
    xdot = y(2);       % Velocidad del carro (m/s)
    alpha = y(3);      % Ángulo del pendulo (rad)
    alphadot = y(4);   % Velocidad angular del pendulo (rad/s)
    
    % Coeficientes de las ecuaciones
    M = Mc + Mp;        % Masa total
    H = Ip + Mp * Lp^2; % Momento de inercia
    F = 0;              % Fuerzas externas
    
    cos_alpha = cos(alpha);
    sin_alpha = sin(alpha);
    
    % Fuerzas del carro y pendulo
    numerator_x = (F + Mp * Lp * sin_alpha * alphadot^2 - Mp * Lp * q * cos_alpha * sin_alpha - Bc * xdot - Bp * Lp * alphadot * cos_alpha);
    denominator_x = (M - Mp^2 * Lp^2 * cos_alpha^2 / H);
    xddot = numerator_x / denominator_x;
    
    numerator_alpha = ((Mp * Lp * q * sin_alpha - Bp * alphadot - Mp * Lp * cos_alpha * xddot));
    alphaddot = numerator_alpha / H;
    
    % Derivadas
    dydt = [xdot;
            xddot;
            alphadot;
            alphaddot];
end
