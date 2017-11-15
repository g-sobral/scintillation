A
B
C
a
b

theta

G = (a*b)/(sqrt((A*C)-(B^2/4))*cos(theta));

phi

V_px
V_py
V_pz

V_D

V_sx = -(V_px-tan(theta)*cos(phi)*V_pz);
V_sy = V_D - (V_py-tan(theta)*sin(phi)*V_pz);

V_eff = sqrt( (C*V_sx^2 - B*V_sx*V_sy + A*V_sy^2) / ((A*C)-(B^2/4)));

L_0

q_0 = 2*pi/L_0;

p

nu = p/2;

Ck_L

Cs_L = Ck_L/((1000/(2*pi))^(2*nu -1));

re
lambda
csi

tmp = besselk((nu-1/2),(q_0*csi))/(2*pi*gamma(nu-1/2));
R_deltaPhi = @(csi) (re^2)*(lambda^2)*sec(theta)*G*Cs_L*(abs(csi/(2*q_0))^(nu-1/2))*tmp;

Z

F2 = lambda*Z*sec(theta)/(2*pi);

f = @(csi, k) 2*R_deltaPhi(csi) - R_deltaPhi(csi-k*F2) - R_deltaPhi(csi+k*F2);
g = @(csi, k) exp(f(csi,k)-f(0,k)) - exp(-f(0,k));

%% I(f)
k = 2*pi*f/Veff;

auxInt = @(csi) g(csi,k)*cos(csi*k);
I = 2*integer(auxInt,0,inf);


