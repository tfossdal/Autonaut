function [xdot, M, foilLimits] = autonaut(x,u,t,V_c,beta_c,V_wind, beta_wind,wave_omega,wave_amp,wave_dir)

% x = [ x y z phi theta psi u v w p q r delta_r fluid_memory(1:4) ϑ1 ϑ2 ϑ3 ϑ1_dot ϑ2_dot ϑ3_dot ]'

if nargin == 0
    x = zeros(29,1); u = zeros(2,1); t=0;
end
if nargin < 4
    V_c = 0; beta_c = 0;
end
if nargin < 6
    V_wind = 0; beta_wind = 0;
end
if nargin < 8
    wave_omega = 0;                 % Wave frequency [rad/s]
    wave_amp = 0;                   % Significant wave height [m]
    wave_dir = deg2rad(0);          % Wave direction [rad]
end

persistent thetas_ddot;  % Persistent variable to solve algebraic loop
if isempty(thetas_ddot)
    thetas_ddot = zeros(3, 1);
end

persistent tau_foil;    % Persistent variable to solve algebraic loop
if isempty(tau_foil)
    tau_foil = zeros(3, 1);
end

% Example of logging data
persistent log;
tf = 200;  % Final time for logging

% Initialize persistent log variable on the first call
if isempty(log)
    log = struct('t', [], 'thetas_ddot', [], 'thetas_dot', [], 'Q_N', [], 'Q_A', [], 'Q_inertia', [], 'tau_foil', [], 'tau_rudder', []);
end


%% Main data
% Environmental data
env.g   = 9.81;                 % acceleration of gravity (m/s^2)
env.rho = 1025;                 % density of water
env.rho_a = 1.225;              % density of air (kg/m^3)

% Vessel data
ANaut.m = 250.0;                  % mass (kg)
ANaut.Loa = 5;                    % Length over all [m]
ANaut.L = 5;                      % Length [m]
ANaut.Lpp = 4.6;                  % Length between perpendiculars [m]
ANaut.B = 1.08;                   % beam (m)
ANaut.T = 0.3;                    % draft (m)
ANaut.CB = 0.505;                       % Block coefficient [-]
ANaut.Cwp = 0.8;                        % Waterplane area coefficient [-]
ANaut.GMT = 0.113;                      % Transverse metacentric height [m]
ANaut.R44 = 0.4 * ANaut.B;              % radii of gyrations (m)
ANaut.R55 = 0.25 * ANaut.L;
ANaut.R66 = 0.25 * ANaut.L;
ANaut.J66 = ANaut.m*ANaut.R66^2;        % Inertia DOF 6 [kg m^2]
ANaut.T_sway = 1;                 % time constant in sway (s)
ANaut.T_yaw = 1;                  % time constant in yaw (s)
ANaut.Umax = 6 * 0.5144;          % 6 knots maximum forward speed (m/s)

% Manoeuvering model
a1.m = ANaut.m;                         % mass [kg]
a1.A11 = 27.7;                          % Added mass DOF 1 [kg]
a1.A22 = 176.9;                         % Added mass DOF 2 [kg]
% a1.A22 = 212.9;                         % Added mass DOF 2 [kg]
a1.A66 = 593.5-ANaut.J66;               % Added inertia DOF 6 [kg m^2]
a1.A26 = 0;                             % Added mass from DOF 6 in DOF 2 [kg m]
a1.A62 = 0;                             % Added mass from DOF 2 in DOF 6 [kg m]
a1.B11 = 138.85;                        % Damping DOF 1 [xx]
a1.B22 = 115.73;                        % Damping DOF 2 [xx]
a1.B66 = 197.8;                         % Damping DOF 6 [xx]
a1.Bv11 = a1.B11/1;                     % Quadratic damping DOF 1
a1.Bv22 = a1.B22/1;                     % Quadratic damping DOF 2
a1.Bv66 = 0;
% a1.Bv11 = 347.13;                     % Quadratic damping DOF 1
% a1.Bv22 = 289.33;                     % Quadratic damping DOF 2
% a1.Bv66 = 0;                          % Quadratic damping DOF 6
% Unsure why these have not been used, but they are in the report

a1.xg = 0;                              % Not measured
a1.J66 = ANaut.J66;                     % System inertia DOF 6 [kg m^2]

% Rudder:
rudd.area = 0.11;                       % Rudder area [m^2]
rudd.ratio = 1.68;                      % Aspect ratio [-]
rudd.xH = -0.4;                         % Lateral force coordinate [m]
rudd.xR = -2.3;                         % Longitudinal rudder position [m]
rudd.CN = 1.56;                         % Rudder coefficient [-]
rudd.tR = 0.3;                          % Drag coefficient [-]
rudd.aH = 0.2;                          % Force factor [-]
rudd.Ts = 0.2;

% Wind model:
wind.Cx = 0.50;                         % Wind coefficient X-direction [-]
wind.Cy = 0.70;                         % Wind coefficient Y-direction [-]
wind.Cn = 0.05;                         % Wind coefficient N-moment [-]
wind.AFw = 0.195;                       % Frontal projected area [m^2]
wind.ALw = 1.5;                         % Lateral projected area [m^2]
wind.Loa = ANaut.Loa;                   % Length over all [m]

% Wave model:
wave_phase = pi/2;                      % Wave phase

% Seakeeping model:
a2.L = 6.146;                           % Length* [m]
a2.B = 1.236;                           % Breadth* [m]
a2.T = 0.287;                           % Draft* [m]
a2.Cwp = 0.683;                         % Waterplane area coefficient [-]
a2.GMT = 0.845;                         % Transverse metacentric height [m]
a2.delta = 0.149;                       % Geometric parameter hull [-]
a2.mu = 0.110;                          % Ratio viscous damping to critical damping [-]
a2.d = 0.1705;                          % Longitudinal distance IMU and CGx [m]

a2.filterGain = 10;                     % Filter gain for phase of pitch response

% IMO Resolution A.685(17)
GM_Lnum = 1/12 * a2.L^2/a2.T;
% Geometric coefficient
Cgeo = 0.373 + 0.023*(a2.B/a2.T);

T4 = 2*a2.B*Cgeo/sqrt(a2.GMT);

% Propulsion model:
a3.mF1 = 3.5;                             % Mass [kg]
a3.mF2 = 1.75;
a3.mF3 = 1.75;
a3.cg = 0.309;                          % Distance from pivot to foil CG [m]
a3.J1 = 0.1091;                          % Calculated from the point that the foils have uniform density rho_f = 500
a3.J2 = a3.J1/2;
a3.J3 = a3.J1/2;
a3.AF1 = 0.250;                         % Area [m^2]
a3.AF2 = 0.125;
a3.AF3 = 0.125;
a3.S1 = 1.3;                            % Span [m]
a3.S2 = 0.65;
a3.S3 = 0.65;
a3.th1max = deg2rad(50);                % Maximum deflection angle [rad]
a3.th2max = deg2rad(45);
a3.th3max = deg2rad(45);
% a3.lambdaF = 5;                         % Wave filter parameter for theta
a3.c_m = 0.192;                         % Mean chord length
a3.c_max = 0.2275;                      % Max chord length
a3.mu = 0.4;                            % Geometric constant for c_tip
a3.c_tip = (a3.c_m - a3.c_max*a3.mu)/(1-a3.mu);
% a3.k_s = 1226;                          % Soft spring constant fore   [N/m]
% a3.k_s = 2200;                          % Medium spring constant fore [N/m]
a3.k_s = 5400;                          % Stiff spring constant fore  [N/m]
a3.offset = 0.0;                        % Offset angle foils [rad]
a3.d = 4.7e-3;                          % Wire diameter torsion spring [m]
a3.D = 6*a3.d;                          % Mean diameter torsion spring [m]
a3.x_s = 60e-3;                         % Distance from pivot point to spring [m]
a3.E = 200e9;                           % Youngs modulus torsion spring [Pa]
a3.N_a = 10.3;                          % Number of active turns for torsion spring
a3.eta = 1;
% a3.zetaBf = 1.2836;                     % Relative damping factor for foil system
% a3.zetaBa = 0.8203;
a3.zetaBf = 3;                          % Relative damping factor for foil system
a3.zetaBa = 3;
a3.alphaS = 12.0;                       % "Stall angle" [deg]
a3.LS = 10;                             % "Stall transition" [deg]
a3.CLs = 0.653;                          % Detached lift coefficient
% a3.Calpha = 0.65;                       % Coefficient for unattached flow affecting circulatory moment
% a3.Cf = 1.0;                            % Coefficient for unsteady flow HC affecting circulatory moment
a3.tF = 0.2;                            % Added resistance/reduction from foils
a3.CDs = 0.286;                          % Detached drag coefficient
a3.pivot = 0.2089;                      % Distance from LE to pivot point [per chord length]
a3.TB = 0.02;                           % Damping time constant foils
a3.TBB = 0.1;

% Spring parameters
spring.E = 200e9;                       % Young's modulus steel
spring.d = 4.7e-3;                      % Spring diameter [m]
spring.D = spring.d*6;                  % Spring mean diameter [m]

% Foil pivot points in body-fixed frame
a3.x1_p = 2.4;                           % Longitudinal position foil 1 [m]
a3.y1_p = 0.0;                           % Lateral position foil 1 [m]
a3.z1_p = 0.7;                           % Vertical position foil 1 [m]
a3.r1_p = [a3.x1_p; a3.y1_p; a3.z1_p];

a3.x2_p = -2.3;                         % Longitudinal position foil 2 [m]
a3.y2_p = -0.3;                         % Lateral position foil 2 [m]
a3.z2_p = 0.7;                          % Vertical position foil 2 [m]
a3.r2_p = [a3.x2_p; a3.y2_p; a3.z2_p];

a3.x3_p = -2.3;                         % Longitudinal position foil 3 [m]
a3.y3_p = 0.3;                          % Lateral position foil 3 [m]
a3.z3_p = 0.7;                          % Vertical position foil 3 [m]
a3.r3_p = [a3.x3_p; a3.y3_p; a3.z3_p];

%% State and current variables

eta = x(1:6);                           % positions
nu = x(7:12);                           % velocity vectors

U = sqrt(nu(1)^2 + nu(2)^2);  % speed
u_c = V_c * cos(beta_c - eta(6));       % current surge velocity
v_c = V_c * sin(beta_c - eta(6));       % current sway velocity
nu_c = [u_c v_c 0 0 0 0 ]';                     % current velocity vector

nu_r = nu - nu_c;                       % relative velocity vector
nu_dot = zeros(6,1);


theta = eta(5);

thetas = x(18:20);                  % foil angles ϑ
thetas_dot = x(21:23);              % foil angular velocities \dot{ϑ}

% Enforce foil angle limits
if thetas(1) >= a3.th1max
    thetas(1) = a3.th1max;
    disp('Beginning')
    if thetas_dot(1) > 0
        disp(['thetas_dot(1): ', num2str(thetas_dot(1))]);
        thetas_dot(1) = 0;
        disp(['thetas_dot(1): ', num2str(thetas_dot(1))]);
    end
    if thetas_ddot(1) > 0
        disp(['thetas_ddot(1): ', num2str(thetas_ddot(1))]);
        thetas_ddot(1) = 0;
        disp(['thetas_ddot(1): ', num2str(thetas_ddot(1))]);
    end
end
if thetas(1) <= -a3.th1max
    thetas(1) = -a3.th1max;
    if thetas_dot(1) < 0
        disp(['thetas_dot(1): ', num2str(thetas_dot(1))]);
        thetas_dot(1) = 0;
        disp(['thetas_dot(1): ', num2str(thetas_dot(1))]);
    end
    if thetas_ddot(1) < 0
        disp(['thetas_ddot(1): ', num2str(thetas_ddot(1))]);
        thetas_ddot(1) = 0;
        disp(['thetas_ddot(1): ', num2str(thetas_ddot(1))]);
    end
end


% if abs(thetas(1)) > a3.th1max
%     thetas(1) = a3.th1max*sign(thetas(1));
%     thetas_dot(1) = 0;
%     thetas_ddot(1) = 0;
% end

%% Manoeuvering subsystem

nuR = [nu_r(1); nu_r(2); nu_r(6)];
u_r = nuR(1); v_r = nuR(2); r = nuR(3);

delta_c = u(1);

% Rudder forces and moments
delta = x(13);  % actual rudder angle
% First order rudder dynamics
delta_dot = (-delta + delta_c) / rudd.Ts;

% Rudder forces
alpha_r = delta - atan2(v_r, u_r);

F_N = 0.5*env.rho*(u_r^2 + v_r^2)*rudd.area*rudd.CN*sin(alpha_r);

tau_rudder = [
    -(1-rudd.tR)*F_N*sin(delta);
    -(1+rudd.aH)*F_N*cos(delta);
    -(rudd.xR + rudd.aH*rudd.xH)*F_N*cos(delta)];
% tau_rudder = [ F_N * (1 - rudd.tR) * sin(delta);
%                F_N * (1 + rudd.aH) * cos(delta);
%                F_N * (rudd.xR + rudd.aH * rudd.xH) * cos(delta) * sin(alpha_r)];

% Wind forces and moments
psi = eta(6);  % heading angle
u_w = V_wind*cos(beta_wind-psi);
v_w = V_wind*sin(beta_wind-psi);

u_rw = nu(1) - u_w;
v_rw = nu(2) - v_w;

gamma_rw = -atan2(v_rw, u_rw);

tau_wind = 0.5*env.rho_a*(u_rw^2 + v_rw^2)*...
    [ - wind.Cx * wind.AFw * cos(gamma_rw);
    wind.Cy * wind.ALw * sin(gamma_rw);
    wind.Cn * wind.ALw * wind.Loa * sin(2*gamma_rw) ];


tau = zeros(3,1) + tau_rudder + tau_wind + tau_foil;
% tau(1) = 200;


% Manoeuvering model
MRB = [ a1.m    0           0;
    0       a1.m        a1.m*a1.xg;
    0       a1.m*a1.xg  a1.J66];
CRB = [ 0               -a1.m*r -a1.m*a1.xg*r;
    a1.m*r          0       0;
    a1.m*a1.xg*r    0       0];

MA = [  a1.A11  0       0;
    0       a1.A22  a1.A26;
    0       a1.A62  a1.A66];
CA = [  0                       0               -a1.A22*v_r-a1.A26*r;
    0                       0               a1.A11*u_r;
    a1.A22*v_r+a1.A26*r     -a1.A11*u_r     0];

% CA = [0 0 0-a1.A26*r;
%       0 0 0;
%       0+a1.A26*r 0 0];

Bp = diag([a1.B11, a1.B22, a1.B66]);

Bv = diag([a1.Bv11*abs(u_r), a1.Bv22*abs(v_r), a1.Bv66*abs(r)]);

M = MRB + MA;
C = CRB + CA;
B = Bp + Bv;            % Blending may be refined, but no linear dampening breaks it

% nu_c_dot = zeros(3,1);              % Current acceleration vector

J = eulerang(eta(4),eta(5),eta(6)); % Kinematic transformation matrix

nuR_dot = M\(tau - C*nuR - B*nuR);
nu_dot(1) = nuR_dot(1);
nu_dot(2) = nuR_dot(2);
nu_dot(6) = nuR_dot(3);
eta_dot = J * nu;

%% Seakeeping subsystem

% RAO model:
C33 = env.rho*env.g*a2.B*a2.L;
C44 = env.rho*env.g*a2.L*a2.B*a2.T*a2.GMT;
C55 = env.rho*env.g*a2.L*a2.B*a2.T*GM_Lnum;

M33 = 2*env.rho*a2.B*a2.T*a2.L;
M44 = (T4/(2*pi))^2*C44;
M55 = 2*env.rho*a2.B*a2.T*a2.L*GM_Lnum;

% Damping approximations
omega = wave_omega;
beta = wave_dir - eta(6);  % wave direction relative to heading
k = omega^2/env.g;
omega_e = abs(omega - k*U*cos(beta));
ke = k*abs(cos(beta));
if omega == 0
    eta_wave = 0;
    R = 0;
    f = 1-k*a2.T;
else
    eta_wave = omega_e/omega;
    R = 2*sin(0.5*k*a2.B*eta_wave^2)*exp(-k*a2.T*eta_wave^2);
    f = sqrt((1-k*a2.T)^2 + (R^2/(k*a2.B*eta_wave^3))^2);
end

kappa = exp(-ke*a2.T);

B33 = a2.L * env.rho*env.g/(k*omega_e) * R^2 / (eta_wave^2);

% ----- B44
% This section if from Tufte (2024), and does not appear in Tufte (2025)

% --- Breadth to draft parameters
Lambda_A = a2.B/a2.T;

B_f = (a2.Cwp-a2.delta)/(1-a2.delta);
% Breadth to draught front
Lambda_Af = B_f/a2.T;

% Calculate the sectional damping
Lambda = [Lambda_A; Lambda_Af];
a = zeros(2,1);
b = zeros(2,1);
d = zeros(2,1);
for i=1:2                           % From Tufte (2024)'s code
    if Lambda(i) >= 1 && Lambda(i) <= 3
        a(i) = -3.94*Lambda(i) + 13.69;
        b(i) = -2.12*Lambda(i) - 1.89;
        d(i) = 1.16*Lambda(i) - 7.97;
    elseif Lambda(i) >= 3 && Lambda(i) <= 6
        a(i) = 0.256*Lambda(i) - 0.286;
        b(i) = -0.11*Lambda(i) - 2.55;
        d(i) = 0.033*Lambda(i) - 1.419;
    else
        fprintf(['Warning: The dimensions for breadth to draught for ' ...
            'the RAO-model is outside scope. Ignoring damping.'])
        a(i) = 0;
        b(i) = 0;
        d(i) = 0;
    end
end

% --- end

b44 = env.rho*R*a2.B^2*sqrt(2*env.g/a2.B)*a.*exp(b.*omega_e^(-1.3)).*omega_e.^(d);

% Integrate sectional damping over length
B44 = a2.L*(a2.delta*b44(1) + (1-a2.delta)*b44(2));

% ----- end

B55 = a2.L*a2.T*GM_Lnum*(env.rho*env.g)/(k*omega_e) * R^2/(eta_wave^2);

Mrao = diag([M33, M44, M55]);
Crao = diag([C33, C44, C55]);
Brao = diag([B33, B44, B55]);

% Wave forces

% Heave force
Z0 = C33*kappa*f*sin(0.5*ke*a2.L)/(0.5*ke*a2.L);
% Z01 = C33*kappa*f*sinc(0.5*ke*a2.L);

% Roll force
K0 = abs(sin(beta))*sqrt(env.rho*env.g^2/omega_e*B44);

% Pitch force
M0 = C33/ke * (sin(0.5*ke*a2.L)/(0.5*ke*a2.L) - cos(0.5*ke*a2.L)) * kappa*f;
% M0 = C33/ke * (sinc(0.5*ke*a2.L) - cos(0.5*ke*a2.L)) * kappa*f;


% Wave elevation
zeta_a = wave_amp;
zeta_Z = zeta_a * cos(omega_e*t);
zeta_M = zeta_a * sin(omega_e*t);
zeta_K = zeta_a * cos(omega_e*t + wave_phase);

tau_wave1 = [Z0*zeta_Z; K0*zeta_K; M0*zeta_M];
if wave_omega == 0
    tau_wave1 = [0;0;0];
    Brao = zeros(size(Brao));
end

% Fluid-memory

% Transversal aspect ratio
Lambda_T = a2.B/a2.T;

% Parameters q0', p0', p1' (Eq. 14a-c)
q0_prime = 0.5696 * (Lambda_T - 0.018) / (Lambda_T + 1.035);
p0_prime = 0.5917 * (Lambda_T - 0.245) / (Lambda_T + 0.612);
p1_prime = 0.7376 * (Lambda_T + 0.394) / (Lambda_T + 0.642);

% This is a wrong one
% Ar = [  0,                          1;
%     -(2*env.g/a2.B) * p1_prime, -(2*env.g/a2.B) * p0_prime  ];
% Br = [0; 1];
% Cr33 = [2*env.rho*env.g*sqrt(2*env.g/a2.B)*(a2.L)*q0_prime, 0];
% Cr55 = [2*env.rho*env.g*sqrt(2*env.g/a2.B)*(a2.L*a2.T*GM_Lnum)*q0_prime, 0];

Ar = [-sqrt(2*env.g/a2.B)*p1_prime, -(2*env.g/a2.B)*p0_prime;
    1,                            0];
Br = [1; 0];
Cr33 = [2*env.rho*env.g*sqrt(2*a2.B*env.g)*(a2.L)*q0_prime, 0];
Cr55 = [2*env.rho*env.g*sqrt(2*a2.B*env.g)*(a2.L*a2.T*GM_Lnum)*q0_prime, 0];



% Dynamics

xi = eta(3:5);
xi_dot = eta_dot(3:5);

xr = x(14:17);   % fluid memory states
xr_dot = zeros(4,1);
% Heave and pitch fluid memory
xr_dot(1:2) = Ar * xr(1:2) + Br * xi_dot(1);
xr_dot(3:4) = Ar * xr(3:4) + Br * xi_dot(3);

mu_r = [Cr33 * xr(1:2);
    0;
    Cr55 * xr(3:4)];

% Seakeeping model
xi_ddot = Mrao\(tau_wave1 - Crao*xi - Brao*xi_dot - mu_r);

phi = eta(4);
theta = eta(5);
% psi = eta(6);


nu_dot(3) = nu_dot(3) + xi_ddot(1)*cos(phi)*cos(theta); % w
nu_dot(4:5) = nu_dot(4:5) + xi_ddot(2:3);           % p,q

%% Propulsion subsystem


thetas_n = thetas + theta;
thetas_n_dot = thetas_dot + nu(5);  % foil angular velocities in NED frame \dot{ϑ_n}
thetas_n_ddot = thetas_ddot + nu_dot(5);  % foil angular accelerations in NED frame \ddot{ϑ_n}


vnb = nu_r(1:3);
vnb_dot = nu_dot(1:3);
wnb = nu_r(4:6);
wnb_dot = nu_dot(4:6);

x_t = x(24:29);  % state vector for Theodorsen
x_t_dot = zeros(6,1);

sigma2 = @(x,x0) 1./(1+exp(-10./deg2rad(5*a3.LS)*(x-x0)));
C3D = @(L) L./(L+2.25);
Ca = @(L) (L-0.26)./(L+0.29);

% Foil #1

Ry1 = [
    cos(thetas(1))  0   sin(thetas(1));
    0               1   0;
    -sin(thetas(1)) 0   cos(thetas(1))];

x__p1 = zeros(3,1);
rpx1 = Ry1 * x__p1;  % position of x__p 1 in body-fixed frame
rbp1 = a3.r1_p;
r_total1 = rbp1 + rpx1;
wbp1 = thetas_dot(1) * [0;1;0];
% wbp_dot1 = thetas_ddot(1) * [0;1;0];                   % PROBLEM, missing theta_ddot
wbp_dot1 = (thetas_ddot(1)+nu_dot(5)) * [0;1;0];                   % PROBLEM, missing theta_ddot

% Foil 1 calculations

% Relative velocity at the foil point, with tan. vel. contributions
vnx1 = vnb + cross(wnb, r_total1) + cross(wbp1, rpx1);

anx1 = vnb_dot + cross(wnb, r_total1) + cross(wnb, cross(wnb, r_total1)) + ...
    2*cross(wnb, cross(wbp1, rpx1)) + cross(wbp_dot1, rpx1) + cross(wbp1, cross(wbp1, rpx1));
if anx1(1) > 1000
    disp(['anx1(1) large: ', num2str(anx1(1))]);
end
alpha1 = atan2(vnx1(3), vnx1(1)) + thetas(1);    % Angle of attack foil 1
% Disp anx1 fpr debugging
% disp(['anx1(3): ', num2str(anx1(3)), ', anx1(1): ', num2str(anx1(1))]);
alpha1_acc = atan2(anx1(3), anx1(1)) + thetas(1);  % Acceleration angle of attack foil 1

U_r1 = sqrt(vnx1(1)^2 + vnx1(3)^2);   % Relative speed foil 1
U_r1_dot = sqrt(anx1(1)^2 + anx1(3)^2); % Relative acceleration foil 1

L1 = a3.S1/a3.c_m;                           % Span foil 1

% Normal moment foil 1

CLn1 = pi*sin(2*alpha1).*(1 - a3.CLs*sigma2(abs(alpha1),deg2rad(a3.alphaS)) + a3.CLs*sigma2(abs(alpha1),pi-deg2rad(a3.alphaS)));
% CDn1 = 2*pi*sin(alpha1).^2.*a3.CDs.*((1-a3.CLs) + ... % 1-a3.CLs might be wrong
CDn1 = 2*pi*sin(alpha1).^2.*a3.CDs.*(...
    a3.CLs.*sigma2(abs(alpha1),deg2rad(a3.alphaS)) - a3.CLs.*sigma2(abs(alpha1),pi-deg2rad(a3.alphaS)));
CNn1 = (CLn1.^2 + CDn1.^2).^0.5*sign(alpha1);

F_N_2Dn1 = 0.5*env.rho*U_r1^2.*a3.c_m*CNn1;

F_N_3Dn1 = F_N_2Dn1*a3.S1*C3D(L1);

x_cp1 = centerOfPressure(alpha1)*a3.c_m;

Q_N_stat1 = F_N_3Dn1 * (x_cp1 - a3.pivot*a3.c_m);

% Foil #2

Ry2 = [
    cos(thetas(2))  0   sin(thetas(2));
    0               1   0;
    -sin(thetas(2)) 0   cos(thetas(2))];

x__p2 = zeros(3,1);
rpx2 = Ry2 * x__p2;  % position of x__p 1 in body-fixed frame
rbp2 = a3.r2_p;
r_total2 = rbp2 + rpx2;
wbp2 = thetas_dot(2) * [0;1;0];
% wbp_dot2 = thetas_ddot(1) * [0;1;0];                   % PROBLEM, missing theta_ddot
wbp_dot2 = (thetas_ddot(1)+nu_dot(5)) * [0;1;0];                   % PROBLEM, missing theta_ddot

% Foil 2 calculations

% Relative velocity at the foil point, with tan. vel. contributions
vnx2 = vnb + cross(wnb, r_total2) + cross(wbp2, rpx2);

anx2 = vnb_dot + cross(wnb, r_total2) + cross(wnb, cross(wnb, r_total2)) + ...
    2*cross(wnb, cross(wbp2, rpx2)) + cross(wbp_dot2, rpx2) + cross(wbp2, cross(wbp2, rpx2));

alpha2 = atan2(vnx2(3), vnx2(1)) + thetas(2);    % Angle of attack foil 2
alpha2_acc = atan2(anx2(3), anx2(1)) + thetas(2);  % Acceleration angle of attack foil 2

U_r2 = sqrt(vnx2(1)^2 + vnx2(3)^2);   % Relative speed foil 2
U_r2_dot = sqrt(anx2(1)^2 + anx2(3)^2); % Relative acceleration foil 2

L2 = a3.S2/a3.c_m;                           % Span foil 2

% Normal moment foil 2

CLn2 = pi*sin(2*alpha2).*(1 - a3.CLs*sigma2(abs(alpha2),deg2rad(a3.alphaS)) + a3.CLs*sigma2(abs(alpha2),pi-deg2rad(a3.alphaS)));
CDn2 = 2*pi*sin(alpha2).^2.*a3.CDs.*((1-a3.CLs) + ... % 1-a3.CLs might be wrong
    a3.CLs.*sigma2(abs(alpha2),deg2rad(a3.alphaS)) - a3.CLs.*sigma2(abs(alpha2),pi-deg2rad(a3.alphaS)));
CNn2 = (CLn2.^2 + CDn2.^2).^0.5*sign(alpha2);

F_N_2Dn2 = 0.5*env.rho*U_r2^2.*a3.c_m*CNn2;

F_N_3Dn2 = F_N_2Dn2*a3.S2*C3D(L2);

x_cp2 = centerOfPressure(alpha2)*a3.c_m;

Q_N_stat2 = F_N_3Dn2 * (x_cp2 - a3.pivot*a3.c_m);

% Foil #1

Ry3 = [
    cos(thetas(3))  0   sin(thetas(3));
    0               1   0;
    -sin(thetas(3)) 0   cos(thetas(3))];

x__p3 = zeros(3,1);
rpx3 = Ry3 * x__p3;  % position of x__p 3 in body-fixed frame
rbp3 = a3.r3_p;
r_total3 = rbp3 + rpx3;
wbp3 = thetas_dot(3) * [0;1;0];
wbp_dot3 = (thetas_ddot(3)+nu_dot(5)) * [0;1;0];                   % PROBLEM, missing theta_ddot

% Foil 3 calculations

% Relative velocity at the foil point, with tan. vel. contributions
vnx3 = vnb + cross(wnb, r_total3) + cross(wbp3, rpx3);

anx3 = vnb_dot + cross(wnb, r_total3) + cross(wnb, cross(wnb, r_total3)) + ...
    2*cross(wnb, cross(wbp3, rpx3)) + cross(wbp_dot3, rpx3) + cross(wbp3, cross(wbp3, rpx3));

alpha3 = atan2(vnx3(3), vnx3(1)) + thetas(3);    % Angle of attack foil 3
alpha3_acc = atan2(anx3(3), anx3(1)) + thetas(3);  % Acceleration angle of attack foil 3

U_r3 = sqrt(vnx3(1)^2 + vnx3(3)^2);   % Relative speed foil 3
U_r3_dot = sqrt(anx3(1)^2 + anx3(3)^2); % Relative acceleration foil 3

L3 = a3.S3/a3.c_m;                           % Span foil 3

% Normal moment foil 3

CLn3 = pi*sin(2*alpha3).*(1 - a3.CLs*sigma2(abs(alpha3),deg2rad(a3.alphaS)) + a3.CLs*sigma2(abs(alpha3),pi-deg2rad(a3.alphaS)));
CDn3 = 2*pi*sin(alpha3).^2.*a3.CDs.*((1-a3.CLs) + ... % 1-a3.CLs might be wrong
    a3.CLs.*sigma2(abs(alpha3),deg2rad(a3.alphaS)) - a3.CLs.*sigma2(abs(alpha3),pi-deg2rad(a3.alphaS)));
CNn3 = (CLn3.^2 + CDn3.^2).^0.5*sign(alpha3);

F_N_2Dn3 = 0.5*env.rho*U_r3^2.*a3.c_m*CNn3;

F_N_3Dn3 = F_N_2Dn3*a3.S3*C3D(L3);

x_cp3 = centerOfPressure(alpha3)*a3.c_m;

Q_N_stat3 = F_N_3Dn3 * (x_cp3 - a3.pivot*a3.c_m);

% Theodorsen model

A_t_11 = [  -1.696*U_r1/a3.c_m, -0.380*(U_r1/a3.c_m)^2;
    1,                  0 ];
B_t_11 = [  1;  0];
C_t_11 = [  0.250*U_r1/a3.c_m,  0.190*(U_r1/a3.c_m)^3 ];
D_t_11 = 0.5;

A_t_22 = [  -1.696*U_r2/a3.c_m, -0.380*(U_r2/a3.c_m)^2;
    1,                  0 ];
B_t_22 = [  1;  0];
C_t_22 = [  0.250*U_r2/a3.c_m,  0.190*(U_r2/a3.c_m)^3 ];
D_t_22 = 0.5;

A_t_33 = [  -1.696*U_r3/a3.c_m, -0.380*(U_r3/a3.c_m)^2;
    1,                  0 ];
B_t_33 = [  1;  0];
C_t_33 = [  0.250*U_r3/a3.c_m,  0.190*(U_r3/a3.c_m)^3 ];
D_t_33 = 0.5;

A_t = blkdiag(A_t_11, A_t_22, A_t_33);
B_t = blkdiag(B_t_11, B_t_22, B_t_33);
C_t = blkdiag(C_t_11, C_t_22, C_t_33);
D_t = blkdiag(D_t_11, D_t_22, D_t_33);

Q_N_stat = [Q_N_stat1; Q_N_stat2; Q_N_stat3];

x_t_dot = A_t*x_t + B_t*Q_N_stat;
Q_N = C_t*x_t + D_t*Q_N_stat;

% Added mass moment
F_A_1 = Ca(L1) * 1/4 * env.rho * pi * a3.c_m^2 * a3.S1 * U_r1_dot * sin(alpha1_acc);
F_A_2 = Ca(L2) * 1/4 * env.rho * pi * a3.c_m^2 * a3.S2 * U_r2_dot * sin(alpha2_acc);
F_A_3 = Ca(L3) * 1/4 * env.rho * pi * a3.c_m^2 * a3.S3 * U_r3_dot * sin(alpha3_acc);
F_A = [F_A_1; F_A_2; F_A_3];
Q_A = F_A*(0.5-0.2089)*a3.c_m;

% Moment inertia
omega = nu(4:6);
omega_dot = nu_dot(4:6);
% Total acceleration of the pivot point in body-fixed frame
rbp1 = a3.r1_p;
r_p1_ddot = nu_dot(1:3) + cross(omega_dot, rbp1) + cross(omega, cross(omega, rbp1));

r_pb1_ddot = Rzyx(eta(4),eta(5),eta(6)) * r_p1_ddot;  % in NED frame
x1_ddot = r_pb1_ddot(1);
z1_ddot = r_pb1_ddot(3);

rbp2 = a3.r2_p;
r_p2_ddot = nu_dot(1:3) + cross(omega_dot, rbp2) + cross(omega, cross(omega, rbp2));

r_pb2_ddot = Rzyx(eta(4),eta(5),eta(6)) * r_p2_ddot;  % in NED frame
x2_ddot = r_pb2_ddot(1);
z2_ddot = r_pb2_ddot(3);

rbp3 = a3.r3_p;
r_p3_ddot = nu_dot(1:3) + cross(omega_dot, rbp3) + cross(omega, cross(omega, rbp3));

r_pb3_ddot = Rzyx(eta(4),eta(5),eta(6)) * r_p3_ddot;  % in NED frame
x3_ddot = r_pb3_ddot(1);
z3_ddot = r_pb3_ddot(3);

% 𝑄 inertia = −𝛿𝑥𝑚 𝑧_𝑝 cos(𝜗𝑛 ) − 𝛿𝑥𝑚𝑥_𝑝 sin(𝜗𝑛 ), <=>
% − Z̈p cos(ϑn ) · (xc.g. − xp )
% − Ẍp sin(ϑn) · (xc.g. − xp )
Q_inertia_1 = (-z1_ddot*cos(thetas_n(1)) -x1_ddot*sin(thetas_n(1)))*a3.mF1* (a3.cg-a3.pivot)*a3.c_m;
Q_inertia_2 = (-z2_ddot*cos(thetas_n(2)) -x2_ddot*sin(thetas_n(2)))*a3.mF2* (a3.cg-a3.pivot)*a3.c_m;
Q_inertia_3 = (-z3_ddot*cos(thetas_n(3)) -x3_ddot*sin(thetas_n(3)))*a3.mF3* (a3.cg-a3.pivot)*a3.c_m;
Q_inertia = [Q_inertia_1; Q_inertia_2; Q_inertia_3];

% Foild dynamics #1

MFF1 = a3.J1 + a3.eta*1/128*env.rho*pi*a3.c_m^4*a3.S1;
MFF2 = a3.J2 + a3.eta*1/128*env.rho*pi*a3.c_m^4*a3.S2;
MFF3 = a3.J3 + a3.eta*1/128*env.rho*pi*a3.c_m^4*a3.S3;
MFF = diag([MFF1, MFF2, MFF3]);

CFF1 = a3.k_s*a3.x_s^2;           % Linear spring; aft foils are different
CFF2 = 1/64*a3.d^4*a3.E/(a3.D*a3.N_a);
CFF3 = 1/63*a3.d^4*a3.E/(a3.D*a3.N_a);
CFF = diag([CFF1, CFF2, CFF3]);

BFFl1 = MFF1/a3.TB;
BFFq1 = MFF1/a3.TBB;
BFF1 = BFFl1 + BFFq1*abs(thetas_n_dot(1));
BFFl2 = MFF2/a3.TB;
BFFq2 = MFF2/a3.TBB;
BFF2 = BFFl2 + BFFq2*abs(thetas_n_dot(2));
BFFl3 = MFF3/a3.TB;
BFFq3 = MFF3/a3.TBB;
BFF3 = BFFl3 + BFFq3*abs(thetas_n_dot(3));
BFF = diag([BFF1, BFF2, BFF3]);

% % Temp
% Q_A_1 = Q_A(1);
% Q_N_1 = Q_N(1);

thetas_n_ddot = MFF\( -BFF*thetas_n_dot - CFF*thetas_n + CFF*ones(3,1)*theta + Q_A + Q_N + Q_inertia );

thetas_ddot = thetas_n_ddot - nu_dot(5);


% Enforce foil angle limits
foilLimits = [a3.th1max, a3.th2max, a3.th3max];
for j = 1:3
    if thetas(j) >= foilLimits(j)
        thetas(j) = foilLimits(j);
        disp('End')
        if thetas_dot(j) > 0
            disp(['thetas_dot(',j,'): ', num2str(thetas_dot(j))]);
            thetas_dot(j) = 0;
            disp(['thetas_dot(',j,'): ', num2str(thetas_dot(j))]);
        end
        if thetas_ddot(j) > 0
            disp(['thetas_ddot(',j,'): ', num2str(thetas_ddot(j))]);
            thetas_ddot(j) = 0;
            disp(['thetas_ddot(',j,'): ', num2str(thetas_ddot(j))]);
        end
    end
    if thetas(j) <= -foilLimits(j)
        thetas(1) = -foilLimits(j);
        if thetas_dot(j) < 0
            disp(['thetas_dot(',j,'): ', num2str(thetas_dot(j))]);
            thetas_dot(j) = 0;
            disp(['thetas_dot(',j,'): ', num2str(thetas_dot(j))]);
        end
        if thetas_ddot(j) < 0
            disp(['thetas_ddot(',j,'): ', num2str(thetas_ddot(j))]);
            thetas_ddot(j) = 0;
            disp(['thetas_ddot(',j,'): ', num2str(thetas_ddot(j))]);
        end
    end
    if abs(thetas_ddot(j)) > 15
        thetas_ddot(j) = 15*sign(thetas_ddot(j));
    end
end

thetas_n = thetas + theta;
F_N = Q_N ./ ([x_cp1;x_cp2;x_cp3] - a3.pivot*a3.c_m);

tau_foil(1) = sum(abs((F_N + F_A) * (1-a3.tF) .* sin(thetas_n)));
if tau_foil(1) > 5000
    disp('Warning: tau_foil is too high')
end

%% Why does theta_n_ddot go to 3000? Debugging
% if thetas_n_ddot(1) > 3000
%     disp('Warning: High foil acceleration detected!');

%     % Components of BFF1
%     disp(['MFF1: ', num2str(MFF1)]);
%     disp(['a3.TB: ', num2str(a3.TB)]);
%     disp(['BFFl1 (MFF1/a3.TB): ', num2str(BFFl1)]);
%     disp(['a3.TBB: ', num2str(a3.TBB)]);
%     disp(['BFFq1 (MFF1/a3.TBB): ', num2str(BFFq1)]);
%     disp(['abs(thetas_dot(1)): ', num2str(abs(thetas_dot(1)))]);
%     disp(['BFF1 (BFFl1 + BFFq1*abs(thetas_dot(1))): ', num2str(BFF1)]);

%     % Components of thetas_ddot(1)
%     disp(['MFF1: ', num2str(MFF1)]);
%     disp(['-BFF1*thetas_dot(1): ', num2str(-BFF1*thetas_dot(1))]);
%     disp(['-CFF1*thetas_n(1): ', num2str(-CFF1*thetas_n(1))]);
%     disp(['CFF1*theta: ', num2str(CFF1*theta)]);
%     disp(['Q_N_1: ', num2str(Q_N_1)]);
%     disp(['Q_A_1: ', num2str(Q_A_1)]);
%     disp(['Q_inertia_1: ', num2str(Q_inertia_1)]);
%     disp(['Sum of components: ', num2str(-BFF1*thetas_dot(1) - CFF1*thetas_n(1) + CFF1*theta + Q_N_1 + Q_A_1 + Q_inertia_1)]);
%     disp(['Resulting thetas_ddot(1): ', num2str(thetas_ddot(1))]);
%     disp('');
% end

% Store thetas_ddot, thetas_dot, and time in the log
log.t(end+1, 1) = t;
log.thetas_dot(end+1, :) = thetas_dot(:)';
log.thetas_ddot(end+1, :) = thetas_ddot(:)';
log.Q_N(end+1, :) = Q_N(:)';
log.Q_A(end+1, :) = Q_A(:)';
log.Q_inertia(end+1, :) = Q_inertia(:)';
log.tau_foil(end+1, 1) = tau_foil(1)';
log.tau_rudder(end+1, :) = tau_rudder(:)';

% Save the log to a MAT-file when the simulation is complete
if t >= tf % Assuming tf is the final simulation time
    save('simulation_log.mat', 'log');
    disp('Simulation log saved to simulation_log.mat');
    disp(t);
end

% % Printing foil accels every integer t for debygging
% if abs(t - round(t)) < 1e-6
%     disp(['Time t: ', num2str(t)]);
%     disp(['thetas_ddot(1): ', num2str(thetas_ddot(1))]);
%     disp(['thetas_ddot(2): ', num2str(thetas_ddot(2))]);
%     disp(['thetas_ddot(3): ', num2str(thetas_ddot(3))]);
%     disp(['tau_foil(1): ', num2str(tau_foil(1))]);
% end



%% Time derivative of the state vector, for numerical integration see SIMautonaut.m
xdot = [eta_dot; nu_dot; delta_dot; xr_dot; thetas_dot; thetas_ddot; x_t_dot];





end