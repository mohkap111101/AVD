function [minValues, minMass] = wingOptimisation(materials, mins, maxs, wing)
    
    options = optimoptions('fmincon','Display','iter','ConstraintTolerance',1e-20,'OptimalityTolerance',1e-20,'UseParallel',true);

    % Setup a global search 

    x_0 = [1.5e-3, 2e-3, 3e-2, 2e-2];

    gs_problem = createOptimProblem('fmincon', ...
        'objective', @(x) wingMass(createGeometry(x(1), x(2), x(3), x(4), wing),  materials), ...
        "x0", x_0, "lb", mins, ....
        "ub", maxs, "options", ...
        options);

    gs = GlobalSearch("Display", "iter", "PlotFcn","gsplotbestf");

    rng(14, "twister")

    [minValues, minMass] = run(gs, gs_problem);
  
end
