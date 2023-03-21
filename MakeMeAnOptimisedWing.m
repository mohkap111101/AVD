% Script to implement a global search to find the optimum values for the
% wing structure 

clear
clc 


load('liftSurfGeom.mat','wing')
material.E = 73.85e9;
material.density = 2.765e3;
material.compYield = 378e6;
wing.L = 0.5;

mins = [0.001, 0.006, 1e-3, 5e-2];
maxs = [10e-2, 10e-2, wing.minHeight, 50e-2];

X_s_0 = [1.5e-2, 1.5e-2, 1e-2, 10e-2];

options = optimoptions('fmincon','Display','iter','ConstraintTolerance',1e-20,'OptimalityTolerance',1e-20,'UseParallel',false);

gs_problem = createOptimProblem('fmincon', ...
        'objective', (@(x) wingMass(createGeometry(x(1), x(2), x(3), x(4), wing), material)), ...
        "x0", X_s_0, "lb", mins, ....
        "ub", maxs, "options", ...
        options);

gs = GlobalSearch("Display","iter", "PlotFcn","gsplotbestf");

rng(14, "twister");

[min_values, min_mass] = run(gs, gs_problem);


rng(14, "twister")

