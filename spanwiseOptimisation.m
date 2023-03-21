% Housekeeping
clear
clc

% Load .mat files
load('liftSurfGeom.mat','wing')
load('rootParams.mat')

% Material properties
density = 2765;
E       = 73.85e9;
sigma   = 378e6;

% Augment wing
wing.b      = b;
wing.h      = h;
wing.t_s    = ts;
wing.A_s    = 1.6 * h * ts;
wing.T      = t + wing.A_s/b;
wing.F      = 0.82;

% Root N
N = realN(wing, 0);
I = inertiaZ(wing);

% Discretise wing
L = sqrt(pi^2 * E * I / (wing.A_s*N/wing.T));

