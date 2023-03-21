function [q02] = Spar_z(geometry, material, z)

% This script is used to determine the web thickness of the wing spars

% N.B. Spars are placed at 20 % and 60 %

% Spar material: Al-7055, T76511
% E=73.5 GPa

E=material.E;
K=8.1;

% The first step is to discretise into stations:

span = geometry.span/2 - geometry.TEfus; 

% Now find the wing torque at each location. 

ds=0.01;
[dist,torque,torque_ult] = RCMK_Spars(ds);

% Now, find the shear flow due to shear force: q2

root_chord=geometry.TEfusC;
tip_chord=geometry.C_r*geometry.taper;
c_r = ((tip_chord-root_chord)/span)*z + root_chord;
coords=geometry.wingboxCoords;

h_r = c_r*((coords(1,2)-coords(2,2))+(coords(4,2)-coords(3,2)))/2;

% Now, find the shear flow due to shear force: q02

dist=round(dist./ds).*ds;

idx = find(dist==round(z/ds)*ds);
torque_ultimate = torque_ult(idx);
q02 = -(torque_ultimate/(2*h_r*0.4*c_r));

end