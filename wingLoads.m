% Housekeeping
clear
clc

% Load lifting surface geometry
load liftSurfGeom.mat

%% Load Parameters
% Load cases
MTOW            = 45119 * 9.81;

%% Wing
% General
wingLoad.disc   = linspace(0,wing.span,100);

% Lift
wingLoad.L_0    = MTOW / (wing.span * 0.672); % L = L_0 * sqrt(1-(y/span)^2)
wingLoad.L      = wingLoad.L_0 * sqrt(1-(wingLoad.disc/wing.span).^2);

% Inertial
