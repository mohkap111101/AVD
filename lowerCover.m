% Housekeeping
clear
clc

% Load wing
load('wingCover.mat')
load('liftSurfGeom.mat','wing')
geometry = wing;

% Material properties
material.density    = 2765;
material.E          = 73.85e9;
material.sigma      = 378e6;

% Initialise Params
Params = [];
goodParams = [];

% Define variables
As_bt   = 0.5;
ts_t    = 0.55;
safety  = 1.15;
sigmaR  = 1.05;

geometry.t      = upper.t;
geometry.b      = upper.b;
geometry.h      = upper.h;
geometry.t_s    = upper.t_s;
geometry.A_s    = 1.6 * upper.h * upper.t_s;
geometry.t_sm   = upper.t + geometry.A_s/upper.b;

% Create lower wing
% Skin thickness distribution
distL = rib.distL;
distT = linspace(0,geometry.span/2-geometry.TEfus);
% distT = distT(distT>2);
% distL = distL(:,distL(1,:)>2);
distT = [distT; LowermodSkin(geometry, material, distT, sigmaR, distL)];
newMass = skinMass(geometry, material, distT);

figure()
plot(distT(1,:),distT(2,:))

%% Save wing
lower.b = geometry.b;
lower.h = geometry.h;
lower.t = geometry.t;
lower.t_s = geometry.t_s;
lower.M = newMass;
lower.distT = distT;
figure()
plot(upper.distT(1,:),upper.distT(2,:))
save('addStraightLines.mat','upper','lower','rib')