% This script analytically calculated the SF and BM 
% distribution about the wing, for Remi

% First, we must define some masses, namely: wing
% structural mass, fuel mass, engine mass, and u/c mass
function BM = Remerz(dist)

wing_mass = 1335.5; % kg
uc_mass = 65; % kg
engine_mass = 3323; % kg 
fuel_mass = 5563.8; % kg

% Next, define where these masses are distributed from:

wing_length = 12.165 - 0; % wing: 0 - 12.165 
uc_length = 0.330 - 0; % u/c: 0 - 0.330
engine_length = 0; % engine: Point at 3.18
fuel_length = 10.95 - 0; % fuel: 0 - 10.95
engine_loc = 3.18;

% Now, find the loads per unit length:
wing_unit = wing_mass / wing_length;
uc_unit = uc_mass / uc_length;
engine_unit = engine_mass;
fuel_unit = fuel_mass / fuel_length;

% Now, start at the end of the wing and find the bending
% moments and shear forces analytically: 

% ----------------------------------
% A   B     C                D     E 

z=wing_length-dist;

% A-B
if dist<uc_length
M = -wing_unit * z * (z/2)...
    - fuel_unit * (z-(wing_length-fuel_length)) * (z-(wing_length-fuel_length))/2 ...
    - engine_unit * (z-(wing_length-engine_loc)) ...
    - uc_unit * (z-(wing_length-uc_length)) * (z-(wing_length-uc_length))/2;
M = M * 9.81;

% B-C
elseif (dist>uc_length) && (dist<engine_loc)
M = -wing_unit * z * (z/2)...
    -fuel_unit * (z-(wing_length-fuel_length)) * (z-(wing_length-fuel_length))/2 ...
    -engine_unit * (z-(wing_length-engine_loc));
M = M * 9.81;

% C-D
elseif (dist>engine_loc) && (dist<fuel_length)
M = -wing_unit * z * (z/2)...
    -fuel_unit* (z-(wing_length-fuel_length)) * (z-(wing_length-fuel_length))/2;
M = M * 9.81;

% D-E
else
M = -wing_unit * z * (z/2);
M = M * 9.81;

end



% Now, output the BM at z:

BM=M;

% Plot the BM distribution: N.B. This is for the left wing:

% figure
% hold on
% fplot(M_DE,[0 wing_length-fuel_length])
% fplot(M_CD,[wing_length-fuel_length wing_length-engine_loc])
% fplot(M_BC,[wing_length-engine_loc wing_length-uc_length])
% fplot(M_AB,[wing_length-uc_length wing_length])
% hold off
end