function [shear, bending] = NEWwingLanding(geometry, y)
% function that takes discretised wing and MTOW and outputs the
% distribution of shear force and bending moments

reaction    = 9.81*2.7*0.85*45119/2;
fuselage    = 3.34/2;
y_gear      = 2-fuselage;

if y<2-geometry.TEfus
    shear = reaction;
    bending = (y-y_gear) * reaction;
else
    shear = 0;
    bending = 0;
end