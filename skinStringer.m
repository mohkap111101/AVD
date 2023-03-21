% Housekeeping
clear
clc

% Load wing geometry
load('liftSurfGeom.mat','wing')
material.E = 73.85e9;
material.density = 2765;
material.compYield = 378e6;
wing.L = 0.5;

% Create arrays for grid search parameter values
h_b     = 0.15:0.05:1; 
b       = 0.05:0.05:0.5;
t_s_t   = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2.0]; % DO NOT CHANGE!!!!!!!
t       = 0.001:0.0005:0.01;

% Initialise bestSoFar
acceptable = 1e5;
goodParams = [];

% Grid search loops
for i = 1:length(h_b)
    for j = 1:length(b)
        for k = 1:length(t_s_t)
            for l = 1:length(t)
                % input parameters to struct
                wing.h      = h_b(i)*b(j);
                wing.b      = b(j);
                wing.t_s    = t_s_t(k)*t(l);
                wing.t      = t(l);

                mass = wingMass(wing, material);
                if mass < acceptable
                    A_s = 1.6 * wing.h * wing.t_s;
                    goodParams = [goodParams; wing.b, wing.h, wing.t,...
                        wing.t_s, h_b(i), t_s_t(k), A_s/wing.b/wing.t, mass];
                end
            end
        end
    end
end

% find the best
[best, index] = min(goodParams(:,7));
goodParams(index,:)