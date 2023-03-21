% This script analytically calculated the SF and BM 
% distribution about the wing, for Remi

% First, we must define some masses, namely: wing
% structural mass, fuel mass, engine mass, and u/c mass
function SF = RoherzHT(dist)

HT_mass = 180; % kg

% Next, define where these masses are distributed from:

HT_length = 9.5380/2 


HT_unit = HT_mass / HT_length;

% Now, start at the end of the wing and find the bending
% moments and shear forces analytically: 

% ----------------------------------
% A   B     C                D     E 

z=HT_length-dist;

% A-B
SF = -HT_unit * z;
SF = SF * 9.81;


% Plot the BM distribution: N.B. This is for the left wing:

% figure
% hold on
% fplot(M_DE,[0 wing_length-fuel_length])
% fplot(M_CD,[wing_length-fuel_length wing_length-engine_loc])
% fplot(M_BC,[wing_length-engine_loc wing_length-uc_length])
% fplot(M_AB,[wing_length-uc_length wing_length])
% hold off
end