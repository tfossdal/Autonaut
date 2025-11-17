function [xdot, M] = autonaut(x,u,t,V_c,beta_c,V_wind, beta_wind,wave_omega,wave_amp,wave_dir)

% x = [ x y z phi theta psi u v w p q r delta_r]' 

if nargin == 0
    x = zeros(17,1); u = zeros(2,1); t=0;
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

T_Nnum = 2*Cgeo*a2.B/sqrt(a2.GMT);

T4 = 2*a2.B*Cgeo/sqrt(a2.GMT);


% State and current variables
eta = x(1:6);                           % positions
nu = x(7:12);                           % velocity vectors

U = sqrt(nu(1)^2 + nu(2)^2);  % speed
u_c = V_c * cos(beta_c - eta(6));       % current surge velocity
v_c = V_c * sin(beta_c - eta(6));       % current sway velocity
nu_c = [u_c v_c 0 0 0 0 ]';                     % current velocity vector

nu_r = nu - nu_c;                       % relative velocity vector
nu_dot = zeros(6,1);

%% Manoeuvering subsystem

nuR = [nu_r(1); nu_r(2); nu_r(6)];
u_r = nuR(1); v_r = nuR(2); r = nuR(3);

delta_c = u(1);

% Rudder forces and moments
delta = x(13);  % actual rudder angle
% First order rudder dynamics
delta_dot = (-delta + delta_c) / rudd.Ts;

% Rudder forces
alpha = delta - atan2(v_r, u_r);

F_N = 0.5*env.rho*(u_r^2 + v_r^2)*rudd.area*rudd.CN*sin(alpha);

tau_rudder = [-(1-rudd.tR)*F_N*sin(delta);
              -(1+rudd.aH)*F_N*cos(delta);
              -(rudd.xR + rudd.aH*rudd.xH)*F_N*cos(delta)];
% tau_rudder = [ F_N * (1 - rudd.tR) * sin(delta);
%                F_N * (1 + rudd.aH) * cos(delta);
%                F_N * (rudd.xR + rudd.aH * rudd.xH) * cos(delta) * sin(alpha)];

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
 

tau = zeros(3,1) + tau_rudder + tau_wind;   
tau(1) = 200;


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
B = Bp + Bv;            % Blending may be refined, but no linear breaks it

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
for i=1:2
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
% Z0 = C33*kappa*f*sinc(0.5*ke*a2.L);

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
q0_prime = 0.5696 * (Lambda_T^2) + 1.035 * (Lambda_T - 0.018);
p0_prime = 0.5917 * (Lambda_T^2) + 0.245 * (Lambda_T - 0.612);
p1_prime = 0.7376 * (Lambda_T^2) + 0.394 * (Lambda_T - 0.642);

Ar = [  0,                          1;
        -(2*env.g/a2.B) * p1_prime, -(2*env.g/a2.B) * p0_prime  ];
Br = [0; 1];
Cr33 = [2*env.rho*env.g*sqrt(2*env.g/a2.B)*(a2.L)*q0_prime, 0];
Cr55 = [2*env.rho*env.g*sqrt(2*env.g/a2.B)*(a2.L*a2.T*GM_Lnum)*q0_prime, 0];



% Dynamics

xi = eta(3:5);
xi_dot = eta_dot(3:5);

xr = x(14:17);   % fluid memory states
xr_dot = zeros(4,1);
% Heave fluid memory
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


nu_dot(3) = nu(3) + xi_ddot(1)*cos(phi)*cos(theta); % w
nu_dot(4:5) = nu_dot(4:5) + xi_ddot(2:3);      % p,q


% Time derivative of the state vector, numerical integration see SIMautonaut.m  
xdot = [eta_dot; nu_dot; delta_dot; xr_dot];





end