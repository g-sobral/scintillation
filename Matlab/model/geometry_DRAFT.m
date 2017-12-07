
% Step 6 - Satellite to xyzp
Rp = [cosd(splong) sind(splong) 0;
    -sind(splong) cosd(splong) 0;
    0 0 1];

uvw = Rp*[xR; yR; zR];

uvw_satp = Rp*[xpL; ypL; zpL];

DC3p = [1 0 0; 
    0 cosd(-splat) sind(-splat);
    0 -sind(-splat) cosd(-splat)];

DCp = DC3p*DC2*DC1;

RLSp = DCp*[uvw(1) - uvw_satp(1)
    uvw(2) - uvw_satp(2)
    uvw(3) - uvw_satp(3)];

xp = RLSp(2);
yp = RLSp(1);
zp = -RLSp(3);

% Calculate Line of Sight
uvw = R*[xS; yS; zS];
RLS = DC*[uvw(1) - xRuvw; uvw(2) - yRuvw; uvw(3) - zRuvw];
xL = RLS(1);
yL = RLS(2);
zL = RLS(3);

deltaT = t(k)-t(k-1); %s

vxL = (xL(k)-xL(k-1))/deltaT;
vyL = (yL(k)-yL(k-1))/deltaT;
vzL = (zL(k)-zL(k-1))/deltaT;

rangeL = sqrt(xL2 + yL^2 + zL^2);

rangep = sqrt(xp^2 + yp^2 + zp^2);

vsf = rangep/rangeL;

V_px = vyL*vsf;

V_py = vxL*vsf;

V_pz = -vzL*vsf;

Theta = rad2deg(atan2(sqrt(xp^2 + yp^2),zp));

Phi = rad2deg(atan2(yp,xp));
