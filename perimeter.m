function coords = perimeter(geometry)
% Takes geometry (struct with lifting surface geometric parameters)
% Outputs x and y values of corners

% Simplify parameters
C_r_LE  = geometry.C_r*(1-(1-geometry.taper)*geometry.LEfus/geometry.span);
C_r_TE  = geometry.C_r*(1-(1-geometry.taper)*geometry.TEfus/geometry.span);
C_t     = geometry.C_r * geometry.taper;
sweep   = geometry.span/2 * sind(geometry.sweep);

% Find coordinates of corners
rootTE  = [geometry.TEfus,  - 0.75 * C_r_TE];
tipTE   = [geometry.span/2, - 0.75 * C_t - sweep];
tipLE   = [geometry.span/2, 0.25 * C_t - sweep];
rootLE  = [geometry.LEfus,  0.25 * C_r_LE];

% Combine into useful array
coords  = [rootTE;
            tipTE;
            tipLE;
            rootLE;
            rootTE];
end