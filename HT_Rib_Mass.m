function [mass] = HT_Rib_Mass (geometry, y, t_eff, spacing)

% Spar material: Al-7055, T76511
% Density = 2850
density = 2850;
E=7.35*10^10;

coords=geometry.wingboxCoords;
local_chord = geometry.C_r*(1-(1-geometry.taper)*(geometry.TEfus + y)/(geometry.span/2));
h = local_chord*((coords(1,2)-coords(2,2))+(coords(4,2)-coords(3,2)))/2;

BM_i=HTweight(geometry,y)*4.5;
BM_a=HTaeroMom(geometry,y);
BM = BM_a - BM_i;
if BM<0
    BM = 0;
end

I = (local_chord * t_eff^3)/12 + (local_chord * t_eff * (h/2)^2);
F = (BM^2 * spacing * h * t_eff * local_chord) / (2*E*I^2);
t_r = ((F * h^2)/(3.62*E*local_chord)).^(1/3);

if t_r<0.001
    t_r=0.001;
end

mass = density * 0.5 * local_chord * h * t_r;