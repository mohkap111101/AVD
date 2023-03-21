function [M , t_E, t_R] = ribThickness(geometry, material, y, BM)

syms x;

c = geometry.C_r*(1-(1-geometry.taper)*(y+geometry.TEfus)/(geometry.span/2));


%% RC's Remerz function to find BM at each rib:
% Bending Moments:

% Now, it is possible to find the wing box height 
fSparHeight = geometry.wingboxCoords(1,2) - geometry.wingboxCoords(2,2);
rSparHeight = geometry.wingboxCoords(4,2) - geometry.wingboxCoords(3,2);
D = (fSparHeight + rSparHeight) / 2;

%% Crushing 

% To determine the rib thickness, it is necessary to equate the yield
% stress of the material to the critical buckling stress of the material

% From Zahra, the critical buckling stress is given by: 3.62E(t/h)^2 where
% E is Young's Modulus, t is rib thickness, and h is wing box depth

% From Daqing, the critical yield stress is given by: F/t*c where F is the
% crush force (see video 11:40), and c is chord at rib position.

% F = (M^2 * s * h * t_e * c)/ (2 * E * I^2) 

syms t_e

t_r = (geometry.L*t_e - geometry.L*geometry.t_sm)/D;

cWingbox=c*geometry.wingboxWidth;
% Crushing force:

I = (cWingbox * t_e^3)/12 + (cWingbox * t_e * (D/2)^2);
F = (BM^2 * geometry.L * D * t_e * cWingbox) / (2*material.E*I^2);
t_r2 = ((F * D^2)/(3.62*material.E*cWingbox))^(1/3);

eq = t_r == t_r2;
t_E = vpasolve(eq, t_e);
t_E = double(t_E);
t_R = (geometry.L*t_E - geometry.L*geometry.t_sm)/D;

if t_R < 0.001
    t_R = 0.001;
end

t_E = geometry.t_sm + (D*t_R)/geometry.L;
M = t_E * cWingbox * material.density;
