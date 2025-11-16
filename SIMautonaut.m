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
V_wind = 2.5;                         % Wind velocity [m/s]
beta_wind = deg2rad(-45);                % Cardinal wind direction [rad]

% Initial states
x = zeros(13,1);                 % x = [xn yn zn phi theta psi u v w p q r]'
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
simdata = zeros(nTimeSteps, 17);    % Preallocate table for simulation data

for i = 1:nTimeSteps

    % delta_c = Kp*e_psi + Kd*e_r + Ki*e_int;
    e_psi = ssa(psi_d - x(6)); % Heading error
    e_psi_int = e_psi_int + e_psi*h; % Integral of heading error
    e_psi_int = sat(e_psi_int, deg2rad(20)/Ki); % Anti-windup for integral term

    delta_c = Kp*(e_psi) + Ki*(e_psi_int) - Kd*x(12); % PID
    delta_c = sat(delta_c, deg2rad(30)); % Rudder saturation
    u(1) = delta_c;


    simdata(i, :) = [x', u', e_psi, e_psi_int];   % Store simulation data
    x = rk4(@autonaut, h, x, u, V_c, beta_c, V_wind, beta_wind); % Integrate using RK4 method

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
delta_r = simdata(:,13);
u    = simdata(:,14:15);
e_psi = simdata(:,16);
e_psi_int = simdata(:,17);

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

end

