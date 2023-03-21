% Housekeeping
clear
clc

% Load wing geometry
load('liftSurfGeom.mat', 'HT');
geometry = HT; %%%%% Change for HT

% Material properties
material.density    = 2765;
material.E          = 73.85e9;
material.sigma      = 378e6;

% Initialise Params
Params = [];
goodParams = [];

% Define variables
As_bt       = 0.5;
ts_t        = 0.55;
N           = HTrealN(geometry, 0); %%%%% CHANGE REMERZ
[tempBM,~]  = HTweight(geometry, 0);
BM          = HTaeroMom(geometry, 0) - 4.5 * tempBM; %%%%% CHANGE REMERZ
safety      = 1.15;

%% Grid search for optimal parameters
t = 0.001:0.0001:0.010;
catchpole = [1;...
            1.3];

for j = 1:size(catchpole,2)

    % Ts/t and sigma_cr/sigma
    ts_t = catchpole(1,j);
    sigmaR = catchpole(2,j);

    for i =1:length(t)
        
        % Find t from equating compressive stress from bending to z stringer
        % buckling stress
        T   = t(i) + As_bt*t(i);
        b   = sqrt(sigmaR*3.62*material.E*t(i)^2 * T / (safety*N));
        ts  = ts_t * t(i);
        h   = As_bt*b*t(i)/(1.6*ts);
    
        % Augment wing
        geometry.t      = t(i);
        geometry.b      = b;
        geometry.h      = h;
        geometry.t_s    = ts;
        geometry.A_s    = 1.6 * h * ts;
        geometry.t_sm   = t(i) + geometry.A_s/b;
    
        % Find rib spacing
        I           = inertiaZ(geometry);
        Lstringer   = sqrt(pi^2 * material.E * I / (geometry.A_s*N/geometry.t_sm));
        Lpanel      = sqrt(3.62*material.E*geometry.t_sm^2/(safety*N/geometry.t_sm));
        geometry.L  = min([Lstringer, Lpanel]);
    
        % Mass
        mass = scuffedMass(geometry, material, 0);
    
        % Check valid inputs
        width = min([geometry.wingboxWidth * geometry.TEfusC, geometry.L]);
        gB = 3.62*material.E*(geometry.t_sm/width)^2 > safety*N/geometry.t_sm;
        sB = eulerBuckle(material.E, inertiaZ(geometry), geometry) > N/geometry.t_sm;
        cY = material.sigma > N/geometry.t_sm;
        mH = geometry.h < geometry.minHeight/2;
        
        % Solutions
        Params = [Params; geometry.b, geometry.h, geometry.t, geometry.t_s, geometry.L, mass,...
                    gB, sB, cY, mH];
        [M , t_E, t_R] = ribThickness(geometry, material, 0, BM);
        if gB && sB && cY && mH
            
            % Rib distribution
            distL = [0;geometry.L];
            spanL = geometry.L;
            Mrib = 0;
            while distL(1,end)+spanL < geometry.span/2 - geometry.TEfus
                distL    = [distL, [distL(1,end)+spanL; spanL]];
                sectN   = abs(HTrealN(geometry, distL(1,end)));
                spanL   = findL(geometry, material, I, sectN);
                if spanL == inf
                    spanL = distL(2,end);
                end
                ribMass = HT_Rib_Mass(geometry, distL(1,end), geometry.t_sm, spanL);
                Mrib    = Mrib + ribMass;
            end
    
            % Skin thickness distribution
            distT = linspace(0,geometry.span/2-geometry.TEfus);
            distT = [distT; HTmodSkin(geometry, material, distT, sigmaR, distL)];
            newMass = skinMass(geometry, material, distT);
            MT = newMass + Mrib;
            
            % Tabulate results
            goodParams = [goodParams; geometry.b, geometry.h, geometry.t, geometry.t_s,...
                            geometry.L, mass, newMass, Mrib, MT];
            [size(goodParams,1), ts]
%             input('Press enter')
        end
        
    end
end

%% Save important stuff
HTskin.b = geometry.b;
HTskin.h = geometry.h;
HTskin.t = geometry.t;
HTskin.t_s = geometry.t_s;
HTskin.M = newMass;
HTskin.distT = distT;

HTrib.M = Mrib;
HTrib.distL = distL;

save('tail.mat', 'HTskin', 'HTrib')
