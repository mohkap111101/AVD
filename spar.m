function [fspar, rspar] = spar(geometry)
% Takes geometry (struct with lifting surface geometric parameters)
% Outputs x and y values of spars

% Simplify parameter names
FRfus   = geometry.LEfus - geometry.fspar*(geometry.LEfus-geometry.TEfus);
REfus   = geometry.LEfus - geometry.rspar*(geometry.LEfus-geometry.TEfus);
C_r_f   = geometry.C_r*(1-(1-geometry.taper)*FRfus/(geometry.span/2));
C_r_r   = geometry.C_r*(1-(1-geometry.taper)*REfus/(geometry.span/2));
C_t     = geometry.C_r * geometry.taper;

% Sweep
sweepFR    = FRfus * sind(geometry.sweep);
sweepRE    = REfus * sind(geometry.sweep);
sweepTip   = geometry.span/2 * sind(geometry.sweep);

% Find coordinates of spar ends
rootF   = [FRfus,           (0.25-geometry.fspar)*C_r_f - sweepFR];
tipF    = [geometry.span/2, (0.25-geometry.fspar)*C_t - sweepTip];
rootR   = [REfus,           (0.25-geometry.rspar)*C_r_r - sweepRE];
tipR    = [geometry.span/2, (0.25-geometry.rspar)*C_t - sweepTip];

% Combine into useful arrays
fspar   = [rootF;
            tipF];

rspar   = [rootR;
            tipR];

end