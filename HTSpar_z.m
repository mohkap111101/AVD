function [q02] = HTSpar_z(geometry, material, z)

% This script is used to determine the web thickness of the wing spars

% N.B. Spars are placed at 20 % and 60 %

% Spar material: Al-7055, T76511
% E=73.5 GPa

E=material.E;
K=8.1;

% The first step is to discretise into stations:

span = geometry.span/2 - geometry.TEfus; 

% Now find the wing torque at each location. 

ht = linspace(0,span,1000);
lift = HTAero(ht, 33000);
w_r = 2 * 360*9.81/((1+geometry.taper)*geometry.span);
weight  = 4.5 * w_r * (1 - (1-geometry.taper) * (geometry.TEfus + ht)/(geometry.span/2));
torque_ult = HTTorque(ht,weight,lift,z);

% Now, find the shear flow due to shear force: q2

root_chord=geometry.TEfusC;
tip_chord=geometry.C_r*geometry.taper;
c_r = ((tip_chord-root_chord)/span)*z + root_chord;
coords=geometry.wingboxCoords;

h_r = c_r*((coords(1,2)-coords(2,2))+(coords(4,2)-coords(3,2)))/2;

% Now, find the shear flow due to shear force: q02

q02 = -(torque_ult/(2*h_r*0.5*c_r));

end