close all;
clear all;
clc;

%% parameters from Carrano 2012
params.Z = 350e3;       % altitude of the phase screen [m]
params.L_0 = 10e3;      % turbulence outer scale [m]
params.a = 50;          % scaling factor a
params.b = 1;           % scaling factor b
params.delta = 0;       % irregularities inclination angle from xz plane

%% constants
params.c = 299792458;   % Speed of Light ~ 300*10^6 [m/s]
params.sf = 1575e6;     % Signal Frequency [Hz]
params.lambda = params.c/params.sf; % radio wavelength [m]
params.r_e = 2.8179403227e-15;  % classical electron radius [m]

%% propagation geometry parameters [TO DO]
params.theta = deg2rad(30);  % propagation (nadir) angle at the IPP
params.phi = deg2rad(30);    % magnetic meridian angle of the propagation vector
params.psi = deg2rad(30);    % magnetic inclination angle
params.V_px = 1;        % IPP velocity component in the magnetic north direction
params.V_py = 1;        % IPP velocity component in the magnetic east direction
params.V_pz = 1;        % IPP velocity component in the down direction

%% independent parameters
params.Ck_L = 5.8e35;   % vert. int. turb strength @ 1km scale [INDEPENDENT]
params.p = 3.39;        % phase spectral index [INDEPENDENT]
params.V_D = 147.4;     % zonal irregularity drift velocity [INDEPENDENT]

[params.V_eff, params.G] = propagation_geo(params);

%%

f_min = 0;
f_max = 25;

f = f_min:0.01:f_max;
I = zeros(size(f));

for i = 1:length(f)
    I(i) = temporal_psd(f(i), params);
end

semilogx(f, 10*log10(I))
