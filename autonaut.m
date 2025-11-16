function [xdot, M] = autonaut(x,u,V_c, beta_c,wind,waves)

% x = [ x y z phi theta psi u v w p q r ]' 

if nargin == 0
    x = zeros(12,1); u = zeros(1,1); V_c = 0; beta_c = 0;
end
    

% Main data
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


% State and current variables
eta = x(1:6);                           % positions
nu = x(7:12);                           % velocity vectors

U = sqrt(nu(1)^2 + nu(2)^2);  % speed
u_c = V_c * cos(beta_c - eta(6));       % current surge velocity
v_c = V_c * sin(beta_c - eta(6));       % current sway velocity
nu_c = [u_c v_c 0 0 0 0 ]';                     % current velocity vector

nu_r = nu - nu_c;                       % relative velocity vector
nu_dot = zeros(6,1);

% Manoeuvering subsystem

nuR = [nu_r(1); nu_r(2); nu_r(6)];
u_r = nuR(1); v_r = nuR(2); r = nuR(3);

delta_c = u(1);

% Rudder forces and moments
delta = delta_c;

alpha = delta - atan2(v_r, u_r);

F_N = 0.5*env.rho*(u_r^2 + v_r^2)*rudd.area*rudd.CN*sin(alpha);

tau_rudder = [-(1-rudd.tR)*F_N*sin(delta);
              -(1+rudd.aH)*F_N*cos(delta);
              -(rudd.xR + rudd.aH*rudd.xH)*F_N*cos(delta)];
% tau_rudder = [ F_N * (1 - rudd.tR) * sin(delta);
%                F_N * (1 + rudd.aH) * cos(delta);
%                F_N * (rudd.xR + rudd.aH * rudd.xH) * cos(delta) * sin(alpha)];

tau = zeros(3,1) + tau_rudder;
tau(1) = 10;


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
B = Bp + Bv;

% nu_c_dot = zeros(3,1);              % Current acceleration vector

J = eulerang(eta(4),eta(5),eta(6)); % Kinematic transformation matrix

nuR_dot = M\(tau - C*nuR - B*nuR);
nu_dot(1) = nuR_dot(1);
nu_dot(2) = nuR_dot(2); 
nu_dot(6) = nuR_dot(3);
eta_dot = J * nu;


% Time derivative of the state vector, numerical integration see SIMautonaut.m  
xdot = [eta_dot; nu_dot];





end