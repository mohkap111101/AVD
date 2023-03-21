function moment = MoKapFailed(z_sub)
% This script analytically calculated the SF and BM 
% distribution about the wing, for Remi

% First, we must define some masses, namely: wing
% structural mass, fuel mass, engine mass, and u/c mass

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

% Locations
A_loc   = wing_length;
B_loc   = wing_length-uc_length;
C_loc   = wing_length-engine_length;
D_loc   = wing_length-fuel_length;
E_loc   = 0;

% Now, start at the end of the wing and find the bending
% moments and shear forces analytically: 

% ----------------------------------
% A   B     C                D     E 

syms z;
% D-E
M_DE = -wing_unit * z * (z/2);
M_DE = M_DE * 9.81;

% C-D
M_CD = -wing_unit * z * (z/2)...
    -fuel_unit* (z-(wing_length-fuel_length)) * (z-(wing_length-fuel_length))/2;
M_CD = M_CD * 9.81;

% B-C
M_BC = -wing_unit * z * (z/2)...
    -fuel_unit * (z-(wing_length-fuel_length)) * (z-(wing_length-fuel_length))/2 ...
    -engine_unit * (z-(wing_length-engine_loc));
M_BC = M_BC * 9.81;

% A-B
M_AB = -wing_unit * z * (z/2)...
    - fuel_unit * (z-(wing_length-fuel_length)) * (z-(wing_length-fuel_length))/2 ...
    - engine_unit * (z-(wing_length-engine_loc)) ...
    - uc_unit * (z-(wing_length-uc_length)) * (z-(wing_length-uc_length))/2;
M_AB = M_AB * 9.81;

% Now, find the maximum BM:

if z_sub < E_loc
    moment = NaN;
elseif z_sub < D_loc
    moment = double(subs(M_DE,z,z_sub));
elseif z_sub < C_loc
    moment = double(subs(M_CD,z,z_sub));
elseif z_sub < B_loc
    moment = double(subs(M_BC,z,z_sub));
elseif z_sub <= A_loc
    moment = double(subs(M_AB,z,z_sub));
else
    moment = NaN;
end

end