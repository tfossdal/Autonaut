function SIMautonaut()
% SIMautonaut - A simple simulation of an autonaut vessel
% This function simulates the dynamics of an autonaut vessel

clearvars;                                  % Clear variables from memory
close all;                                  % Close all figure windows
clear autonaut                              % Clear persistent variables

%% USER INPUTS
h  = 0.01;                       % Sampling time [s]
T_final = 100;	                 % Final simulation time [s]

psi_d_deg = 10;                     % Desired heading (deg)
psi_d = deg2rad(psi_d_deg);         % Desired heading (rad)

% Ocean current
V_c = 0.3*0;                       % Ocean current speed (m/s)
beta_c = deg2rad(0);            % Ocean current direction (rad)

% Initial states
x = zeros(12,1);                 % x = [xn yn zn phi theta psi u v w p q r]'
u = zeros(2,1);                 % Control input vector, u = [ delta_c thrust_c]

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

[~, M] = autonaut();    % Get mass matrix M

% PID heading autopilot parameters (Nomoto model: M(6,6) = T/K)
T = 1;                           % Nomoto time constant
% K = T / M(6,6);                  % Nomoto gain constant
K = T / M(3,3);                  % Nomoto gain constant

wn = 1.5;                        % Closed-loop natural frequency (rad/s)
zeta = 1.0;                      % Closed-loop relative damping factor (-)

% Kp = M(6,6) * wn^2;                     % Proportional gain
Kp = M(3,3) * wn^2;                     % Proportional gain
% Kd = M(6,6) * (2 * zeta * wn - 1/T);    % Derivative gain
Kd = M(3,3) * (2 * zeta * wn - 1/T);    % Derivative gain
Td = Kd / Kp;                           % Derivative time constant
Ti = 10 / wn;                           % Integral time constant

%% MAIN LOOP
simdata = zeros(nTimeSteps, 14);    % Preallocate table for simulation data

for i = 1:nTimeSteps

    % delta_c = Kp*e_psi + Kd*e_r + Ki*e_int;
    delta_c = Kp*(ssa(psi_d - x(6))) - Kd*x(12); % No integral action
    delta_c = sat(delta_c, deg2rad(30)); % Rudder saturation
    u(1) = delta_c;


    simdata(i, :) = [x', u'];   % Store simulation data
    x = rk4(@autonaut, h, x, u, V_c, beta_c); % Integrate using RK4 method

    % --- Progress Update Logic ---
    if mod(i, floor(nTimeSteps / 10)) == 0 % Print an update every 10%
        progress = (i/nTimeSteps) * 100;
        fprintf('  %d%% complete\n', round(progress));
    elseif i == 1
        disp("  Simulation in progress:")
        fprintf('  %d%% complete\n', 0);
    end
end

%% PLOTS
eta  = simdata(:,1:6); 
nu   = simdata(:,7:12); 
u    = simdata(:,13:14);

figure()
plot(eta(:,2), eta(:,1), 'linewidth', 2);
legend('Autonaut Path');
title('Autonaut Path in NED frame');
xlabel('Y (m)');
ylabel('X (m)');

figure()
plot(t, [psi_d_deg*ones(nTimeSteps,1), rad2deg(eta(:,6))], 'linewidth', 2);
legend('Desired Heading', 'Actual Heading');
title('Autonaut Heading');
xlabel('Time (s)');
ylabel('Heading (deg)');

figure()
plot(t, rad2deg(u(:,1)), 'linewidth', 2);
title('Rudder Command');
xlabel('Time (s)');
ylabel('Rudder Angle (deg)');

