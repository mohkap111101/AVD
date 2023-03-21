function [shear, bending] = wingLanding(wing, MTOW)
% function that takes discretised wing and MTOW and outputs the
% distribution of shear force and bending moments

reaction    = 9.81*2.7*0.85*MTOW/2;
fuselage    = 3.34/2;
y_gear      = 2-fuselage;

shear               = zeros(size(wing));
shear(wing<=y_gear) = reaction;

bending                 = zeros(size(wing));
bending(wing<y_gear)    = (y_gear-wing(wing<y_gear))*reaction;