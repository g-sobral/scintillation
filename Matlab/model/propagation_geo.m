function [V_eff, G] = propagation_geo(params)
a = params.a;           % scaling factor
b = params.b;           % scaling factor
theta = params.theta;   % propagation (nadir) angle at the IPP
phi = params.phi;       % magnetic meridian angle of the propagation vector
psi = params.psi;       % magnetic inclination angle
delta = params.delta;   % irregularities inclination angle from xz plane
V_px = params.V_px;     % IPP velocity component in the magnetic north direction
V_py = params.V_py;     % IPP velocity component in the magnetic east direction
V_pz = params.V_pz;     % IPP velocity component in the down direction
V_D = params.V_D;       % zonal irregularity drift velocity [INDEPENDENT]

%% trigonometric constants
c_theta = cos(theta);
t_theta = tan(theta);
s_phi = sin(phi);
c_phi = cos(phi);
s_psi = sin(psi);
c_psi = cos(psi);
s_delta = sin(delta);
c_delta = cos(delta);

%% propagation coefficients
C11 = (a^2)*(c_psi^2) + (s_psi^2)*((b^2)*(s_delta^2) + (c_delta^2));
C22 = (b^2)*(c_delta^2) + (s_delta^2);
C33 = (a^2)*(s_psi^2) + (c_psi^2)*((b^2)*(s_delta^2) + (c_delta^2));
C12 = ((b^2) -1)*s_psi*s_delta*c_delta;
C13 = ((a^2) - (b^2)*(s_delta^2) - (c_delta^2))*s_psi*c_psi;
C23 = -((b^2) - 1)*c_psi*s_delta*c_delta;

A = C11 + C33*(t_theta^2)*(c_phi^2) - 2*C13*t_theta*c_phi;
B = 2*(C12 + C33*(t_theta^2)*s_phi*c_phi - t_theta*(C13*s_phi + C23*c_phi));
C = C22 + C33*(t_theta^2)*(s_phi^2) - 2*C23*t_theta*s_phi;

%% geometric enhancement factor
G = (a*b) / (sqrt((A*C) - (B^2)/4)*c_theta);

%% effective scan velocity
V_sx = -(V_px - (t_theta*c_phi*V_pz));
V_sy = V_D - (V_py - (t_theta*s_phi*V_pz));

V_eff = sqrt( (C*(V_sx^2) - B*V_sx*V_sy + A*(V_sy^2)) / (A*C - (B^2)/4));
