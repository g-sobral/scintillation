function [I] = temporal_psd(f, params)
% TEMPORAL_PSD  Temporal power spectral densities. (Carrano 2012,
%               equations 11 to 14)
%

Z       = params.Z;         % altitude of the phase screen [m]
lambda  = params.lambda;    % radio wavelength [m]
theta   = params.theta;     % propagation (nadir) angle at the IPP
V_eff   = params.V_eff;     % effective scan velocity [m/s]

k = 2*pi*f/V_eff;           % spatial wavenumber

F2 = lambda*Z*sec(theta)/(2*pi);    % squared Fresnel parameter

% auxiliary functions, depending on csi
% where csi is the separation distance in the transverse plane
f = @(csi) 2*r_delta_phi((csi), params) - r_delta_phi((csi-k*F2), params) -...
    r_delta_phi((csi+k*F2), params);
g = @(csi) exp(-f(0)-f(csi)) - exp(-f(0));
h = @(csi) g(csi).*cos(csi*k);

% numerically integrates function h(csi) from 0 to Inf
I = 2*integral(h,0,Inf);

function [R] = r_delta_phi(csi, params)
% R_DELTA_PHI  Correlation function of phase fluctuations just
%              beneath the screen. (Carrano 2012, equation 2)
%

% [TO DO] besselk function returns Inf for csi=0
%         besselk function returns i for csi<0  
csi = abs(csi);
if csi == 0
    R = Inf;
    return
end

p       = params.p;         % phase spectral index [INDEPENDENT]
G       = params.G;         % geometric enhancement factor
Ck_L    = params.Ck_L;      % vert. int. turb strength @ 1km scale [INDEPENDENT]
L_0     = params.L_0;       % turbulence outer scale [km]
lambda  = params.lambda;    % radio wavelength
theta   = params.theta;     % propagation (nadir) angle at the IPP
r_e     = params.r_e;       % classical electron radius

q_0 = 2*pi/L_0;             % outer scale wavenumber
nu = p/2;                   % irregularity spectral index

% turbulence strength (Cs) * layer thickness (L)
Cs_L = Ck_L / ((1000/(2*pi))^(2*nu -1));

aux1 = G*Cs_L*(abs(csi/(2*q_0)).^(nu-1/2));
aux2 = besselk((nu-1/2),(q_0.*csi))/(2*pi*gamma(nu-1/2));
R = (r_e^2)*(lambda^2)*sec(theta).*aux1.*aux2;