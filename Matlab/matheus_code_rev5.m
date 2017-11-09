clear
close all
clc

% Input & Fixed Parameters

for runType = 1:6
    
    if runType == 1
        FORT = 1;
        LRange = 1;
        prefix = 'FORT\ResultsL1';
    elseif runType == 2
        FORT = 1;
        LRange = 2;
        prefix = 'FORT\ResultsL2';
    elseif runType == 3
        FORT = 1;
        LRange = 5;
        prefix = 'FORT\ResultsL5';
    elseif runType == 4
        FORT = 0;
        LRange = 1;
        prefix = 'POAL\ResultsL1';
    elseif runType == 5
        FORT = 0;
        LRange = 2;
        prefix = 'POAL\ResultsL2';
    else
        FORT = 0;
        LRange = 5;
        prefix = 'POAL\ResultsL5';
    end
    
    if FORT
        ReceiverLat = -3.7446;           % Station Latitude [deg]
        ReceiverLong = -38.5776;         % Station Longitude [deg]
        ReceiverH = 22.2489;             % Station Altitude [m]
    else
        ReceiverLat = -30.0739;          % Station Latitude [deg]
        ReceiverLong = -51.1197;         % Station Longitude [deg]
        ReceiverH = 77.1771;             % Station Altitude [m]
    end
    
    H = 350e3;                           % Height of Scattering Layer [m]
    
    if LRange == 1
        f = 1575e6;                      % Signal Frequency [Hz]
    elseif LRange == 2
        f = 1227e6;                      % Signal Frequency [Hz]
    else
        f = 1176e6;                      % Signal Frequency [Hz]
    end
    
    f0 = 0;                              % Outer Scale Frequency [Hz]
    a = 50;                              % Scaling Factor
    b = 1;                               % Scaling Factor
    delta = 0;                           % Irregularities Inclin from xz [rad]
    c = 299792458;                       % Speed of Light ~ 300*10^6 m/s
    kfs = 2*pi*f/c;                      % Free Space Signal Wavenumber [1/m]
    tc = 10;                             % Time constant of the HP filter [sec]
    
    % Initialization, Conversion Rates and MATLAB Parameters
    wgs = wgs84Ellipsoid('meters');
    Re = earthRadius('meters');
    deg2semi = 1/180;
    semi2deg = 180;
    
    % Supplied Indexes
    dayLTIdx = 2;
    satidIdx = 5;
    
    if LRange == 1
        s4Idx = 10;
        sigPhiIdx = 16;
        pIdx = 33;
    elseif  LRange == 2
        s4Idx = 35;
        sigPhiIdx = 41;
        pIdx = 47;
    else
        s4Idx = 49;
        sigPhiIdx = 55;
        pIdx = 61;
    end
    
    
    hourIdx = 67;
    minIdx = 68;
    psiIdx = 75;
    xSIdx = 92;
    ySIdx = 93;
    zSIdx = 94;
    
    % Calculated Indexes
    xLIdx = 96;
    yLIdx = 97;
    zLIdx = 98;
    rangeLIdx = 99;
    uxLIdx = 100;
    uyLIdx = 101;
    uzLIdx = 102;
    xpLIdx = 103;
    ypLIdx = 104;
    zpLIdx = 105;
    splatIdx = 106;
    splongIdx = 107;
    elevCIdx = 108;
    azimCIdx = 109;
    earthCAIdx = 110;
    ipplatCIdx = 111;
    ipplongCIdx = 112;
    xpIdx = 113;
    ypIdx = 114;
    zpIdx = 115;
    rangepIdx = 116;
    uxpIdx = 117;
    uypIdx = 118;
    uzpIdx = 119;
    vxLIdx = 120;
    vyLIdx = 121;
    vzLIdx = 122;
    vmagLIdx = 123;
    vxLIdx = 124;
    vyLIdx = 125;
    vzLIdx = 126;
    vmagLIdx = 127;
    vsfIdx = 128;
    vpxIdx = 129;
    vpyIdx = 130;
    vpzIdx = 131;
    thetapCIdx = 132;
    phipCIdx = 133;
    AIdx = 134;
    BIdx = 135;
    CIdx = 136;
    peIdx = 137;
    GIdx = 138;
    Vsx0Idx = 139;
    Vsy0Idx = 140;
    fresnIdx = 141;
    VeffIdx = 142;
    FspIdx = 143;
    FstIdx = 144;
    Vd0Idx = 145;
    baskIdx = 146;
    VdIdx = 147;
    Vd2Idx = 148;
    
    QIdx = 149;
    VeffInfIdx = 150;
    Vd0InfIdx = 151;
    Vd1InfIdx = 152;
    VdInfIdx = 153;
    
    %% Receiver Related Calculations
    
    [xR,yR,zR] = geodetic2ecef(wgs,ReceiverLat,ReceiverLong, ReceiverH);
    
    R = [cosd(ReceiverLong) sind(ReceiverLong) 0;
        -sind(ReceiverLong) cosd(ReceiverLong) 0;
        0 0 1];
    
    origin_uvw = R*[xR; yR; zR];
    xRuvw = origin_uvw(1);
    yRuvw = origin_uvw(2);
    zRuvw = origin_uvw(3);
    
    DC1 = [1 0 0;0 0 1;0 -1 0];
    
    DC2 = [0 0 -1;0 1 0;1 0 0];
    
    DC3 = [1 0 0;
        0 cosd(-ReceiverLat) sind(-ReceiverLat);
        0 -sind(-ReceiverLat) cosd(-ReceiverLat)];
    
    DC = DC3*DC2*DC1;
    
    if FORT
        files(1).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Fortaleza\gpsFORTsp13n.mat';
        files(1).fileName = 'gpsFORTsp13n.mat';
        files(2).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Fortaleza\gpsFORTsp14n.mat';
        files(2).fileName = 'gpsFORTsp14n.mat';
        files(3).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Fortaleza\gpsFORTsp15n.mat';
        files(3).fileName = 'gpsFORTsp15n.mat';
        files(4).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Fortaleza\gpsFORTsu14n.mat';
        files(4).fileName = 'gpsFORTsu14n.mat';
        files(5).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Fortaleza\gpsFORTsu15n.mat';
        files(5).fileName = 'gpsFORTsu15n.mat';
        files(6).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Fortaleza\gpsFORTsu16n.mat';
        files(6).fileName = 'gpsFORTsu16n.mat';
    else
        files(1).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Porto Alegre\gpsPOALsp13n.mat';
        files(1).fileName = 'gpsPOALsp13n.mat';
        files(2).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Porto Alegre\gpsPOALsp14n.mat';
        files(2).fileName = 'gpsPOALsp14n.mat';
        files(3).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Porto Alegre\gpsPOALsp15n.mat';
        files(3).fileName = 'gpsPOALsp15n.mat';
        files(4).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Porto Alegre\gpsPOALsu13n.mat';
        files(4).fileName = 'gpsPOALsu13n.mat';
        files(5).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Porto Alegre\gpsPOALsu14n.mat';
        files(5).fileName = 'gpsPOALsu14n.mat';
        files(6).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Porto Alegre\gpsPOALsu15n.mat';
        files(6).fileName = 'gpsPOALsu15n.mat';
        files(7).url = 'C:\Users\Matheus\Documents\Mestrado\Scripts\Dados\Porto Alegre\gpsPOALsu16n.mat';
        files(7).fileName = 'gpsPOALsu16n.mat';
    end
    
    
    
    for ii = 1:length(files)
        load(files(ii).url)
        if FORT
            newx = x;
            clear x
        end
        
        %% Cleaning Up Data and Creating Structure
        
        s4LL = 0.3;
        s4HL = 1e5;
        sigphiLL = 0.05;
        sigphiHL = 1e5;
        elevLL = 1e5;
        pLL = 1;
        pHL = 4;
        
        % Cleaning Up Data
        clear x auxData stl_days GPS_Data
        x = newx;
        xClean = x(x(:,s4Idx)> s4LL & x(:,s4Idx)<s4HL,:);
        xClean = xClean(xClean(:,sigPhiIdx)> sigphiLL  & xClean(:,sigPhiIdx)< sigphiHL,:);
        xClean = xClean(xClean(:,8)< elevLL,:);
        xClean = xClean(xClean(:,pIdx)> pLL  & xClean(:,pIdx)< pHL,:);
        
        stl_id = unique(xClean(:,satidIdx));
        
        for i = 1:length(stl_id)
            auxData = xClean(xClean(:,satidIdx)==stl_id(i),:);
            stl_days = unique(auxData(:,dayLTIdx));
            if length(stl_days) == 1
                continue
            else
                jReal = 0;
                for j = 1:length(stl_days)
                    if size(auxData(auxData(:,dayLTIdx)==stl_days(j),:),1) == 1
                        jReal = jReal + 1;
                        continue
                    else
                        GPS_Data(i).date(j-jReal).id = stl_id(i);
                        GPS_Data(i).date(j-jReal).day = stl_days(j);
                        GPS_Data(i).date(j-jReal).data = ...
                            auxData(auxData(:,dayLTIdx)==stl_days(j),:);
                    end
                end
            end
        end
        
        clear auxBool
        for i = 1:length(GPS_Data)
            auxBool(i,1) = ~isempty(GPS_Data(i).date);
        end
        
        GPS_Data = GPS_Data(auxBool);
        
        %% Receiver, Satellite and IPP Geometry Procedure
        
        for i = 1:length(GPS_Data)
            for j = 1:length(GPS_Data(i).date)
                for k = 1:size(GPS_Data(i).date(j).data,1)
                    
                    % Calculate Satellite ECEF
                    %             [GPS_Data(i).date(j).data(k,xSIdx),...
                    %                 GPS_Data(i).date(j).data(k,ySIdx),...
                    %                 GPS_Data(i).date(j).data(k,zSIdx)] =...
                    %                 geodetic2ecef(wgs,GPS_Data(i).date(j).data(k,GPSlatIdx),...
                    %                 GPS_Data(i).date(j).data(k,GPSlongIdx),...
                    %                 GPS_Data(i).date(j).data(k,GPSHIdx));
                    
                    % Calculate Line of Sight
                    uvw = R*[GPS_Data(i).date(j).data(k,xSIdx)
                        GPS_Data(i).date(j).data(k,ySIdx)
                        GPS_Data(i).date(j).data(k,zSIdx)];
                    RLS = DC*[uvw(1) - xRuvw; uvw(2) - yRuvw; uvw(3) - zRuvw];
                    GPS_Data(i).date(j).data(k,xLIdx) = RLS(1);
                    GPS_Data(i).date(j).data(k,yLIdx) = RLS(2);
                    GPS_Data(i).date(j).data(k,zLIdx) = RLS(3);
                    
                    % Calculate Range
                    GPS_Data(i).date(j).data(k,rangeLIdx) = ...
                        sqrt(GPS_Data(i).date(j).data(k,xLIdx)^2 + ...
                        GPS_Data(i).date(j).data(k,yLIdx)^2 + ...
                        GPS_Data(i).date(j).data(k,zLIdx)^2);
                    
                    % Unity Vector
                    GPS_Data(i).date(j).data(k,uxLIdx) = ...
                        GPS_Data(i).date(j).data(k,xLIdx)/ ...
                        GPS_Data(i).date(j).data(k,rangeLIdx);
                    
                    GPS_Data(i).date(j).data(k,uyLIdx) = ...
                        GPS_Data(i).date(j).data(k,yLIdx)/ ...
                        GPS_Data(i).date(j).data(k,rangeLIdx);
                    
                    GPS_Data(i).date(j).data(k,uzLIdx) = ...
                        GPS_Data(i).date(j).data(k,zLIdx)/ ...
                        GPS_Data(i).date(j).data(k,rangeLIdx);
                    
                    % Interception Point Iteration
                    err = 10000;
                    h_test = 0;
                    rng_test_min=0;
                    rng_test_max = GPS_Data(i).date(j).data(k,rangeLIdx);
                    
                    while abs(h_test-H)>err
                        rng_test=(rng_test_min+rng_test_max)/2;
                        
                        tcsx_test = GPS_Data(i).date(j).data(k,uxLIdx)*rng_test;
                        tcsy_test = GPS_Data(i).date(j).data(k,uyLIdx)*rng_test;
                        tcsz_test = GPS_Data(i).date(j).data(k,uzLIdx)*rng_test;
                        
                        rls_test = inv(DC)*[tcsx_test; tcsy_test; tcsz_test];
                        
                        xrls = rls_test(1) + xRuvw;
                        yrls = rls_test(2) + yRuvw;
                        zrls = rls_test(3) + zRuvw;
                        
                        rls_test2 = inv(R)*[xrls; yrls; zrls];
                        
                        xrls = rls_test2(1);
                        yrls = rls_test2(2);
                        zrls = rls_test2(3);
                        
                        [lat,lon,h_test] = ecef2geodetic(wgs,xrls,yrls,zrls);
                        
                        %                 fprintf([num2str(h_test) '\n'])
                        %                 pause
                        
                        if h_test>=H;
                            rng_test_max=rng_test;
                        else
                            rng_test_min=rng_test;
                        end
                    end
                    GPS_Data(i).date(j).data(k,xpLIdx) = xrls;
                    GPS_Data(i).date(j).data(k,ypLIdx) = yrls;
                    GPS_Data(i).date(j).data(k,zpLIdx) = zrls;
                    GPS_Data(i).date(j).data(k,splatIdx) = lat;
                    GPS_Data(i).date(j).data(k,splongIdx) = lon;
                    
                    % Step 3 - Determine Elevation and Azimuth of the Satellite
                    
                    GPS_Data(i).date(j).data(k,elevCIdx) = ...
                        atan2(GPS_Data(i).date(j).data(k,zLIdx), ...
                        sqrt(GPS_Data(i).date(j).data(k,xLIdx)^2 + ...
                        GPS_Data(i).date(j).data(k,yLIdx)^2));
                    GPS_Data(i).date(j).data(k,elevCIdx) = rad2deg(...
                        GPS_Data(i).date(j).data(k,elevCIdx));
                    
                    GPS_Data(i).date(j).data(k,azimCIdx) = ...
                        atan2(GPS_Data(i).date(j).data(k,xLIdx), ...
                        GPS_Data(i).date(j).data(k,yLIdx));
                    GPS_Data(i).date(j).data(k,azimCIdx) = rad2deg(...
                        GPS_Data(i).date(j).data(k,azimCIdx));
                    
                    % Step 4 - Earth Centered Angle
                    
                    GPS_Data(i).date(j).data(k,earthCAIdx) = ...
                        0.0137/(GPS_Data(i).date(j).data(k,elevCIdx)*deg2semi...
                        + 0.11) - 0.022;
                    
                    % Step 5 - IPP Latitude and Longitude
                    
                    GPS_Data(i).date(j).data(k,ipplatCIdx) = ReceiverLat* ...
                        deg2semi + GPS_Data(i).date(j).data(k,earthCAIdx)* ...
                        cosd(GPS_Data(i).date(j).data(k,azimCIdx));
                    if (GPS_Data(i).date(j).data(k,ipplatCIdx) > 0.416)
                        GPS_Data(i).date(j).data(k,ipplatCIdx) = 0.416;
                    elseif (GPS_Data(i).date(j).data(k,ipplatCIdx) < -0.416)
                        GPS_Data(i).date(j).data(k,ipplatCIdx) = -0.416;
                    end
                    GPS_Data(i).date(j).data(k,ipplatCIdx) = semi2deg* ...
                        GPS_Data(i).date(j).data(k,ipplatCIdx);
                    
                    GPS_Data(i).date(j).data(k,ipplongCIdx) = ReceiverLong* ...
                        deg2semi + GPS_Data(i).date(j).data(k,earthCAIdx)* ...
                        sind(GPS_Data(i).date(j).data(k,azimCIdx))/ ...
                        cosd(GPS_Data(i).date(j).data(k,ipplatCIdx));
                    GPS_Data(i).date(j).data(k,ipplongCIdx) = semi2deg* ...
                        GPS_Data(i).date(j).data(k,ipplongCIdx);
                    
                    % Step 6 - Satellite to xyzp
                    
                    Rp = [cosd(GPS_Data(i).date(j).data(k,splongIdx)) ...
                        sind(GPS_Data(i).date(j).data(k,splongIdx)) 0;
                        -sind(GPS_Data(i).date(j).data(k,splongIdx)) ...
                        cosd(GPS_Data(i).date(j).data(k,splongIdx)) 0;
                        0 0 1];
                    
                    uvw = Rp*[xR; yR; zR];
                    
                    uvw_satp = Rp*[GPS_Data(i).date(j).data(k,xpLIdx)
                        GPS_Data(i).date(j).data(k,ypLIdx)
                        GPS_Data(i).date(j).data(k,zpLIdx)];
                    
                    DC3p = [1 0 0; 0 ...
                        cosd(-GPS_Data(i).date(j).data(k,splatIdx)) ...
                        sind(-GPS_Data(i).date(j).data(k,splatIdx)); 0 ...
                        -sind(-GPS_Data(i).date(j).data(k,splatIdx)) ...
                        cosd(-GPS_Data(i).date(j).data(k,splatIdx))];
                    
                    DCp = DC3p*DC2*DC1;
                    
                    RLSp = DCp*[uvw(1) - uvw_satp(1)
                        uvw(2) - uvw_satp(2)
                        uvw(3) - uvw_satp(3)];
                    
                    GPS_Data(i).date(j).data(k,xpIdx) = RLSp(2);
                    GPS_Data(i).date(j).data(k,ypIdx) = RLSp(1);
                    GPS_Data(i).date(j).data(k,zpIdx) = -RLSp(3);
                    
                    GPS_Data(i).date(j).data(k,rangepIdx) = ...
                        sqrt(GPS_Data(i).date(j).data(k,xpIdx)^2 + ...
                        GPS_Data(i).date(j).data(k,ypIdx)^2 + ...
                        GPS_Data(i).date(j).data(k,zpIdx)^2);
                    
                    GPS_Data(i).date(j).data(k,uxpIdx) = ...
                        GPS_Data(i).date(j).data(k,xpIdx)/ ...
                        GPS_Data(i).date(j).data(k,rangepIdx);
                    GPS_Data(i).date(j).data(k,uypIdx) = ...
                        GPS_Data(i).date(j).data(k,ypIdx)/ ...
                        GPS_Data(i).date(j).data(k,rangepIdx);
                    GPS_Data(i).date(j).data(k,uzpIdx) = ...
                        GPS_Data(i).date(j).data(k,zpIdx)/ ...
                        GPS_Data(i).date(j).data(k,rangepIdx);
                    
                end
            end
        end
        
        % Step 3.1 - Local Velocity
        
        for i = 1:length(GPS_Data)
            for j = 1:length(GPS_Data(i).date)
                for k = 2:size(GPS_Data(i).date(j).data,1)
                    deltaT = ((GPS_Data(i).date(j).data(k,hourIdx)-...
                        GPS_Data(i).date(j).data(k-1,hourIdx))*60 +...
                        GPS_Data(i).date(j).data(k,minIdx)-...
                        GPS_Data(i).date(j).data(k-1,minIdx))*60;
                    GPS_Data(i).date(j).data(k,vxLIdx) =...
                        (GPS_Data(i).date(j).data(k,xLIdx)-...
                        GPS_Data(i).date(j).data(k-1,xLIdx))/deltaT;
                    GPS_Data(i).date(j).data(k,vyLIdx) =...
                        (GPS_Data(i).date(j).data(k,yLIdx)-...
                        GPS_Data(i).date(j).data(k-1,yLIdx))/deltaT;
                    GPS_Data(i).date(j).data(k,vzLIdx) =...
                        (GPS_Data(i).date(j).data(k,zLIdx)-...
                        GPS_Data(i).date(j).data(k-1,zLIdx))/deltaT;
                    GPS_Data(i).date(j).data(k,vmagLIdx) =...
                        sqrt(GPS_Data(i).date(j).data(k,vxLIdx)^2 +...
                        GPS_Data(i).date(j).data(k,vyLIdx)^2 +...
                        GPS_Data(i).date(j).data(k,vzLIdx)^2);
                end
            end
        end
        
        
        
        for i = 1:length(GPS_Data)
            for j = 1:length(GPS_Data(i).date)
                for k = 1:size(GPS_Data(i).date(j).data,1)
                    
                    % Step 3.5 - IPP Velocity - Vsat_TCS
                    
                    GPS_Data(i).date(j).data(k,vsfIdx) = ...
                        GPS_Data(i).date(j).data(k,rangepIdx)/ ...
                        GPS_Data(i).date(j).data(k,rangeLIdx);
                    
                    GPS_Data(i).date(j).data(k,vpxIdx) = ...
                        GPS_Data(i).date(j).data(k,vyLIdx)* ...
                        GPS_Data(i).date(j).data(k,vsfIdx);
                    
                    GPS_Data(i).date(j).data(k,vpyIdx) = ...
                        GPS_Data(i).date(j).data(k,vxLIdx)* ...
                        GPS_Data(i).date(j).data(k,vsfIdx);
                    
                    GPS_Data(i).date(j).data(k,vpzIdx) = ...
                        -GPS_Data(i).date(j).data(k,vzLIdx)* ...
                        GPS_Data(i).date(j).data(k,vsfIdx);
                    
                    GPS_Data(i).date(j).data(k,thetapCIdx) = rad2deg( ...
                        atan2(sqrt(GPS_Data(i).date(j).data(k,xpIdx)^2 + ...
                        GPS_Data(i).date(j).data(k,ypIdx)^2), ...
                        GPS_Data(i).date(j).data(k,zpIdx)));
                    
                    GPS_Data(i).date(j).data(k,phipCIdx) = rad2deg( ...
                        atan2(GPS_Data(i).date(j).data(k,ypIdx), ...
                        GPS_Data(i).date(j).data(k,xpIdx)));
                    
                    % Geometry Parameters (A, B, C, pe and G)
                    
                    C11 = a^2*cosd(GPS_Data(i).date(j).data(k,psiIdx))^2+...
                        sind(GPS_Data(i).date(j).data(k,psiIdx))...
                        ^2*(b^2*sin(delta)^2+cos(delta)^2);
                    C22 = b^2*cos(delta)^2+sin(delta)^2;
                    C33 = a^2*sind(GPS_Data(i).date(j).data(k,psiIdx))^2+...
                        cosd(GPS_Data(i).date(j).data(k,psiIdx))^2*(b^2*...
                        sind(GPS_Data(i).date(j).data(k,psiIdx))^2 + cos(delta)^2);
                    C12 = (b^2-1)*sind(GPS_Data(i).date(j).data(k,psiIdx))*...
                        sin(delta)*cos(delta);
                    C13 = (a^2-b^2*sin(delta)^2-cos(delta)^2)*...
                        sind(GPS_Data(i).date(j).data(k,psiIdx))*...
                        cosd(GPS_Data(i).date(j).data(k,psiIdx));
                    C23 = -(b^2-1)*cosd(GPS_Data(i).date(j).data(k,psiIdx))*...
                        sin(delta)*cos(delta);
                    
                    GPS_Data(i).date(j).data(k,AIdx) = C11 + C33*...
                        tand(GPS_Data(i).date(j).data(k,thetapCIdx))^2*...
                        cosd(GPS_Data(i).date(j).data(k,phipCIdx))^2 - 2*C13*...
                        tand(GPS_Data(i).date(j).data(k,thetapCIdx))*...
                        cosd(GPS_Data(i).date(j).data(k,phipCIdx));
                    
                    GPS_Data(i).date(j).data(k,BIdx) = 2*(C12 + C33*...
                        tand(GPS_Data(i).date(j).data(k,thetapCIdx))^2*...
                        sind(GPS_Data(i).date(j).data(k,phipCIdx))*...
                        cosd(GPS_Data(i).date(j).data(k,phipCIdx))...
                        -tand(GPS_Data(i).date(j).data(k,thetapCIdx))*(C13*...
                        sind(GPS_Data(i).date(j).data(k,phipCIdx))+C23*...
                        cosd(GPS_Data(i).date(j).data(k,phipCIdx))));
                    
                    GPS_Data(i).date(j).data(k,CIdx) = C22 + C33*...
                        tand(GPS_Data(i).date(j).data(k,thetapCIdx))^2*...
                        sind(GPS_Data(i).date(j).data(k,phipCIdx))^2 - 2*C23*...
                        tand(GPS_Data(i).date(j).data(k,thetapCIdx))*...
                        sind(GPS_Data(i).date(j).data(k,phipCIdx));
                    
                    A_l = (GPS_Data(i).date(j).data(k,AIdx)*...
                        cosd(GPS_Data(i).date(j).data(k,phipCIdx))^2+...
                        GPS_Data(i).date(j).data(k,BIdx)*...
                        cosd(GPS_Data(i).date(j).data(k,phipCIdx))*...
                        sind(GPS_Data(i).date(j).data(k,phipCIdx))+...
                        GPS_Data(i).date(j).data(k,CIdx)*...
                        sind(GPS_Data(i).date(j).data(k,phipCIdx))^2)*...
                        cosd(GPS_Data(i).date(j).data(k,thetapCIdx))^2;
                    
                    B_l = (GPS_Data(i).date(j).data(k,BIdx)*...
                        cosd(2*GPS_Data(i).date(j).data(k,phipCIdx))+...
                        (GPS_Data(i).date(j).data(k,CIdx)-...
                        GPS_Data(i).date(j).data(k,AIdx))*...
                        sind(2*GPS_Data(i).date(j).data(k,phipCIdx)))*...
                        cosd(GPS_Data(i).date(j).data(k,thetapCIdx));
                    
                    C_l = GPS_Data(i).date(j).data(k,AIdx)*...
                        sind(GPS_Data(i).date(j).data(k,phipCIdx))^2 - ...
                        GPS_Data(i).date(j).data(k,BIdx)*...
                        cosd(GPS_Data(i).date(j).data(k,phipCIdx))*...
                        sind(GPS_Data(i).date(j).data(k,phipCIdx)) + ...
                        GPS_Data(i).date(j).data(k,CIdx)*...
                        cosd(GPS_Data(i).date(j).data(k,phipCIdx))^2;
                    
                    D_l = sqrt((C_l-A_l)^2 + B_l^2);
                    A_dl = (A_l+ C_l + D_l)/2;
                    C_dl = (A_l + C_l - D_l)/2;
                    
                    F_hgeom = hypergeom([(1-GPS_Data(i).date(j).data(k,pIdx))/2,...
                        1/2],1,(A_dl - C_dl)/A_dl);
                    
                    GPS_Data(i).date(j).data(k,peIdx) = (a*b)/(sqrt(A_dl)*...
                        C_dl^(GPS_Data(i).date(j).data(k,pIdx)/2))*F_hgeom;
                    
                    GPS_Data(i).date(j).data(k,GIdx) = (a*b)/ ...
                        (cosd(GPS_Data(i).date(j).data(k,thetapCIdx))* ...
                        sqrt(GPS_Data(i).date(j).data(k,AIdx)* ...
                        GPS_Data(i).date(j).data(k,CIdx) - ...
                        GPS_Data(i).date(j).data(k,BIdx)^2/4));
                    
                    % Step 3.6 - Vsx0 e Vsy0
                    
                    
                    GPS_Data(i).date(j).data(k,Vsx0Idx) =...
                        GPS_Data(i).date(j).data(k,vpxIdx)-...
                        tand(GPS_Data(i).date(j).data(k,thetapCIdx))*...
                        cosd(GPS_Data(i).date(j).data(k,phipCIdx))*...
                        GPS_Data(i).date(j).data(k,vpzIdx);
                    
                    GPS_Data(i).date(j).data(k,Vsy0Idx) =...
                        GPS_Data(i).date(j).data(k,vpyIdx)-...
                        tand(GPS_Data(i).date(j).data(k,thetapCIdx))*...
                        sind(GPS_Data(i).date(j).data(k,phipCIdx))*...
                        GPS_Data(i).date(j).data(k,vpzIdx);
                    
                    % Step 3.7 - Veff
                    
                    GPS_Data(i).date(j).data(k,fresnIdx) =...
                        sqrt(H*secd(GPS_Data(i).date(j).data(k,thetapCIdx))/kfs);
                    
                    GPS_Data(i).date(j).data(k,FspIdx) =...
                        (gamma(5-GPS_Data(i).date(j).data(k,pIdx))/4)/...
                        ((2^((GPS_Data(i).date(j).data(k,pIdx)-1)/2))*...
                        sqrt(pi)*(GPS_Data(i).date(j).data(k,pIdx)-1)*...
                        gamma(GPS_Data(i).date(j).data(k,pIdx)+1)/4);
                    
                    GPS_Data(i).date(j).data(k,FstIdx) =...
                        (sqrt(pi)*gamma(GPS_Data(i).date(j).data(k,pIdx)/2))/...
                        (((2*pi)^((GPS_Data(i).date(j).data(k,pIdx)+1)))*...
                        gamma((GPS_Data(i).date(j).data(k,pIdx)+1)/2));
                    
                    GPS_Data(i).date(j).data(k,Vd0Idx) =...
                        GPS_Data(i).date(j).data(k,Vsy0Idx)-...
                        GPS_Data(i).date(j).data(k,BIdx)/...
                        (2*GPS_Data(i).date(j).data(k,AIdx))*...
                        GPS_Data(i).date(j).data(k,Vsx0Idx);
                    
                    GPS_Data(i).date(j).data(k,VeffIdx) =...
                        (GPS_Data(i).date(j).data(k,fresnIdx)/tc)*...
                        (((GPS_Data(i).date(j).data(k,pIdx)-1)/2)*...
                        (GPS_Data(i).date(j).data(k,FspIdx)/...
                        GPS_Data(i).date(j).data(k,FstIdx))*...
                        (GPS_Data(i).date(j).data(k,peIdx)/...
                        GPS_Data(i).date(j).data(k,GIdx))*...
                        (GPS_Data(i).date(j).data(k,sigPhiIdx)^2/...
                        GPS_Data(i).date(j).data(k,s4Idx)^2))...
                        ^(1/(GPS_Data(i).date(j).data(k,pIdx)-1));
                    
                    GPS_Data(i).date(j).data(k,baskIdx) =...
                        1/GPS_Data(i).date(j).data(k,AIdx)*...
                        sqrt((GPS_Data(i).date(j).data(k,AIdx)*...
                        GPS_Data(i).date(j).data(k,CIdx)-...
                        GPS_Data(i).date(j).data(k,BIdx)^2/4)*...
                        (GPS_Data(i).date(j).data(k,AIdx)*...
                        GPS_Data(i).date(j).data(k,VeffIdx)^2-...
                        GPS_Data(i).date(j).data(k,Vsx0Idx)^2));
                    
                    GPS_Data(i).date(j).data(k,VdIdx) =...
                        GPS_Data(i).date(j).data(k,Vd0Idx)+...
                        GPS_Data(i).date(j).data(k,baskIdx);
                    
                    GPS_Data(i).date(j).data(k,Vd2Idx) =...
                        GPS_Data(i).date(j).data(k,Vd0Idx)-...
                        GPS_Data(i).date(j).data(k,baskIdx);
                    
                    % Simplifed Method
                    
                     GPS_Data(i).date(j).data(k,QIdx) = ...
                         ((GPS_Data(i).date(j).data(k,pIdx)-1)/2*...
                         GPS_Data(i).date(j).data(k,FspIdx)/...
                         GPS_Data(i).date(j).data(k,FstIdx)*...
                         gamma(GPS_Data(i).date(j).data(k,pIdx)/2)/...
                         (sqrt(pi)*...
                         gamma((GPS_Data(i).date(j).data(k,pIdx)+1)/2)))...
                         ^(1/(GPS_Data(i).date(j).data(k,pIdx)-1));
                    
                    GPS_Data(i).date(j).data(k,VeffInfIdx) =...
                        GPS_Data(i).date(j).data(k,fresnIdx)/tc*...
                        GPS_Data(i).date(j).data(k,QIdx)*...
                        ((GPS_Data(i).date(j).data(k,sigPhiIdx)/...
                        GPS_Data(i).date(j).data(k,s4Idx)))...
                        ^(2/(GPS_Data(i).date(j).data(k,pIdx)-1));
                     
                     GPS_Data(i).date(j).data(k,Vd0InfIdx) = ...
                         GPS_Data(i).date(j).data(k,vpyIdx) +...
                         (GPS_Data(i).date(j).data(k,vpxIdx)*...
                         sind(GPS_Data(i).date(j).data(k,psiIdx))-...
                         GPS_Data(i).date(j).data(k,vpzIdx)*...
                         cosd(GPS_Data(i).date(j).data(k,psiIdx)))*...
                         sind(GPS_Data(i).date(j).data(k,phipCIdx))*...
                         tand(GPS_Data(i).date(j).data(k,thetapCIdx))/...
                         (cosd(GPS_Data(i).date(j).data(k,psiIdx))-...
                         cosd(GPS_Data(i).date(j).data(k,phipCIdx))*...
                         sind(GPS_Data(i).date(j).data(k,psiIdx))*...
                         tand(GPS_Data(i).date(j).data(k,thetapCIdx)));
                     
                     GPS_Data(i).date(j).data(k,Vd1InfIdx) = ...
                         sqrt(1+...
                         (sind(GPS_Data(i).date(j).data(k,phipCIdx))^2*...
                         tand(GPS_Data(i).date(j).data(k,thetapCIdx))^2)/...
                         (cosd(GPS_Data(i).date(j).data(k,psiIdx))-...
                         cosd(GPS_Data(i).date(j).data(k,phipCIdx))*...
                         sind(GPS_Data(i).date(j).data(k,psiIdx))*...
                         tand(GPS_Data(i).date(j).data(k,thetapCIdx)))^2)*...
                         GPS_Data(i).date(j).data(k,VeffInfIdx);
                    
                    GPS_Data(i).date(j).data(k,VdInfIdx) = ...
                        GPS_Data(i).date(j).data(k,Vd0InfIdx)+...
                        GPS_Data(i).date(j).data(k,Vd1InfIdx);
                
                end
            end
        end
        
        
        save([prefix '\results' files(ii).fileName],'GPS_Data')
        
    end
end
