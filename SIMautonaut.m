function SIMautonaut()
% SIMautonaut - A simple simulation of an autonaut vessel
% This function simulates the dynamics of an autonaut vessel

clearvars;                                  % Clear variables from memory
close all;                                  % Close all figure windows
clear autonaut                              % Clear persistent variables

%% USER INPUTS
h  = 0.01;                       % Sampling time [s]
T_final = 200;	                 % Final simulation time [s]

psi_d_deg = 10;                     % Desired heading (deg)
psi_d = deg2rad(psi_d_deg);         % Desired heading (rad)

% Ocean current
V_c = 0.3*0;                       % Ocean current speed (m/s)
beta_c = deg2rad(0);            % Ocean current direction (rad)

% Wind
V_wind = 0;                         % Wind velocity [m/s]
beta_wind = deg2rad(-45);                % Cardinal wind direction [rad]

% Waves
wave_omega = 1.5;                     % Wave frequency [rad/s]
wave_height = 1;           % Significant wave height [m]
wave_amp = wave_height/2;           % Wave amplitude [m]
wave_dir = deg2rad(40);              % Wave direction [rad]


% Initial states
x = zeros(29,1);                % x = [xn yn zn phi theta psi u v w p q r delta_r ...
                                %  fluid_memory(4) ϑ(3) ϑ_dot(3) x_t(6) ]'
u = zeros(2,1);                 % Control input vector, u = [ delta_c thrust_c]

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final
nTimeSteps = length(t);         % Number of time steps

[~, ~, foilLimits] = autonaut();    % Get mass matrix M amd foil Limits

% % PID heading autopilot parameters (Nomoto model: M(6,6) = T/K)
% T = 1;                           % Nomoto time constant
% % K = T / M(6,6);                  % Nomoto gain constant
% K = T / M(3,3);                  % Nomoto gain constant

% wn = 1.5;                        % Closed-loop natural frequency (rad/s)
% zeta = 1.0;                      % Closed-loop relative damping factor (-)

% Kp = M(6,6) * wn^2;                     % Proportional gain
% Kp = M(3,3) * wn^2;                     % Proportional gain
% Kd = M(6,6) * (2 * zeta * wn - 1/T);    % Derivative gain
% Kd = M(3,3) * (2 * zeta * wn - 1/T);    % Derivative gain
% Td = Kd / Kp;                           % Derivative time constant
% Ti = 10 / wn;                           % Integral time constant
Kp = 1;
Ki = 0.05;
Kd = 10;
e_psi_int = 0;

%% MAIN LOOP
simdata_x = zeros(nTimeSteps, 29);    % Preallocate table for simulation data
simdata_input = zeros(nTimeSteps, 4);    % Preallocate table for simulation data

for i = 1:nTimeSteps

    % delta_c = Kp*e_psi + Kd*e_r + Ki*e_int;
    e_psi = ssa(psi_d - x(6)); % Heading error
    e_psi_int = e_psi_int + e_psi*h; % Integral of heading error
    e_psi_int = sat(e_psi_int, deg2rad(20)/Ki); % Anti-windup for integral term

    delta_c = Kp*(e_psi) + Ki*(e_psi_int) - Kd*x(12); % PID
    delta_c = sat(delta_c, deg2rad(30)); % Rudder saturation
    u(1) = delta_c;


    simdata_x(i, :) = x';   % Store simulation data
    simdata_input(i, :) = [u', e_psi, e_psi_int];   % Store input data data
    x = rk4(@autonaut, h, x, u, t(i), V_c, beta_c, V_wind, beta_wind, wave_omega, wave_amp, wave_dir); % Integrate using RK4 method

    % Enforcement of bounds
    thetas = x(18:20);                  % foil angles ϑ
    thetas_dot = x(21:23);              % foil angular velocities \dot{ϑ}

    % Enforce foil angle limits
    for j = 1:3
        if thetas(j) >= foilLimits(j)
            thetas(j) = foilLimits(j);
            if thetas_dot(j) > 0
                thetas_dot(j) = 0;
            end
        end
        if thetas(j) <= -foilLimits(j)
            thetas(j) = -foilLimits(j);
            if thetas_dot(j) < 0
                thetas_dot(j) = 0;
            end
        end
    end
    x(18:20) = thetas;                  % foil angles ϑ
    x(21:23) = thetas_dot;              % foil angular velocities \dot{ϑ}


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
eta  = simdata_x(:,1:6);
theta = eta(:,5);               % pitch angle
nu   = simdata_x(:,7:12);
q = nu(:,5);                  % pitch rate
delta_r = simdata_x(:,13);
xr = simdata_x(:,14:17);
thetas = simdata_x(:,18:20);
thetas_dot = simdata_x(:,21:23);
u    = simdata_input(:,1:2);
e_psi = simdata_input(:,3);
e_psi_int = simdata_input(:,4);

figure()
plot(eta(:,2), eta(:,1), 'linewidth', 2);
axis equal;
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
subplot(3,1,1)
plot(t, [rad2deg(u(:,1)), rad2deg(delta_r)], 'linewidth', 2);
title('Rudder Angle');
legend('Commanded Rudder Angle', 'Actual Rudder Angle');
xlabel('Time (s)');
ylabel('Rudder Angle (deg)');

subplot(3,1,2)
plot(t, [rad2deg(e_psi), rad2deg(e_psi_int*Ki)], 'linewidth', 2);
title('Heading Error and Integral of Heading Error');
legend('Heading Error', 'Integral of Heading Error');
xlabel('Time (s)');
ylabel('Error (deg)');

figure()
plot(t, nu(:,1), 'linewidth', 2);
title('Surge Velocity');
xlabel('Time (s)');
ylabel('Surge Velocity (m/s)');

figure()
subplot(2,1,1)
plot(t, [eta(:,3)], 'linewidth', 2);
title('Heave, Roll and Pitch');
xlabel('Time (s)');
ylabel('Heave (m)');
legend('Heave');
subplot(2,1,2)
plot(t, rad2deg(eta(:,4:5)), 'linewidth', 2);
xlabel('Time (s)');
ylabel('Roll and Pitch (deg)');
legend('Roll', 'Pitch');

figure()
plot(t, [xr(:,1), xr(:,3)], 'linewidth', 2);
title('Fluid Memory States');
xlabel('Time (s)');
ylabel('Fluid Memory States');

figure()
subplot(2,1,1)
plot(t, [thetas, theta], 'linewidth', 2);
title('Foil Angles');
legend('Foil 1', 'Foil 2', 'Foil 3', 'Pitch Angle');
xlabel('Time (s)');
ylabel('Foil Angles (rad)');

subplot(2,1,2)
plot(t, [thetas_dot, q], 'linewidth', 2);
title('Foil Angular Velocities');
legend('Foil 1', 'Foil 2', 'Foil 3', 'Pitch Rate');
xlabel('Time (s)');

sim_log = load("simulation_log.mat");
sim_log = sim_log.log;

figure()
plot(sim_log.t, sim_log.thetas_dot, 'linewidth', 2);
title('Logged Foil Angular Velocities');
legend('Foil 1', 'Foil 2', 'Foil 3');
xlabel('Time (s)');
ylabel('Foil Angular Velocities (rad/s)');

figure()
plot(sim_log.t, sim_log.thetas_ddot, 'linewidth', 2);
title('Logged Foil Angular Accelerations');
legend('Foil 1', 'Foil 2', 'Foil 3');
xlabel('Time (s)');
ylabel('Foil Angular Accelerations (rad/s^2)');

figure()
subplot(3,1,1)
plot(sim_log.t, sim_log.Q_N, 'linewidth', 2);
title('Logged Foil Normal Forces');
legend('Foil 1', 'Foil 2', 'Foil 3');
xlabel('Time (s)');
ylabel('Force (m)');

subplot(3,1,2)
plot(sim_log.t, sim_log.Q_A, 'linewidth', 2);
title('Logged Foil Added Mass Forces');
legend('Foil 1', 'Foil 2', 'Foil 3');
xlabel('Time (s)');
ylabel('Force (N)');

subplot(3,1,3)
plot(sim_log.t, sim_log.Q_inertia, 'linewidth', 2);
title('Logged Foil Inertia Forces');
legend('Foil 1', 'Foil 2', 'Foil 3');
xlabel('Time (s)');
ylabel('Moment (Nm)');

figure()
subplot(3,1,1)
plot(sim_log.t, sim_log.tau_foil, 'linewidth', 2);
title('Logged Foil Torques');
xlabel('Time (s)');
ylabel('Torque (Nm)');

subplot(3,1,2)
plot(sim_log.t, sim_log.tau_rudder, 'linewidth', 2);
title('Logged Rudder Forces');
xlabel('Time (s)');
ylabel('Force (N)');
legend('X-direction', 'Y-direction', 'Yaw-direction');



end

