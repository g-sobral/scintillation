function [I] = temporal_psd(f, params)
% TEMPORAL_PSD  Temporal power spectral densities. (Carrano 2012,
%               equations 11 to 14)
%

Z       = params.Z;         % altitude of the phase screen
lambda  = params.lambda;    % radio wavelength
theta   = params.theta;     % propagation (nadir) angle at the IPP

k = 2*pi*f/V_eff;           % spatial wavenumber

F2 = lambda*Z*sec(theta)/(2*pi);    % squared Fresnel parameter

% auxiliary functions, depending on csi
% where csi is the separation distance in the transverse plane
f = @(csi) 2*r_delta_phi(csi) - r_delta_phi(csi-k*F2) - r_delta_phi(csi+k*F2);
g = @(csi) exp(f(csi)-f(0)) - exp(-f(0));
h = @(csi) g(csi)*cos(csi*k);

% numerically integrates function h(csi) from 0 to Inf
I = 2*integral(h,0,Inf);

function [R] = r_delta_phi(csi, params)
% R_DELTA_PHI  Correlation function of phase fluctuations just
%              beneath the screen. (Carrano 2012, equation 2)
%

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

aux1 = G*Cs_L*(abs(csi/(2*q_0))^(nu-1/2));
aux2 = besselk((nu-1/2),(q_0*csi))/(2*pi*gamma(nu-1/2)); 
R = (r_e^2)*(lambda^2)*sec(theta)*aux1*aux2;