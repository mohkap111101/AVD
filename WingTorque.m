function [comb_torque, torque] = WingTorque(lift_per_span, wing_weight_per_span, fuel_weight_per_span, engine_weight_load, uc_weight_per_span, chord_dist, y, x_ac, M_0_torque, Thrust, c_bar_wing, ctip, c0, b)
% Function to determine the torque distributions and total torque along the
% span of the wing
% Outputs -> y              array of stations along the wing span 
%            comb_torque    array of the torques at each station along the
%                           span
%            torque         array of the total torque at each station
%                           along the span

load Mohkap_workspace.mat;

ds = y(2)-y(1);
y = round(y./ds)*ds;

%% Define position of shear centre with respect to chord and the fuselage nose

shear_centre_over_c = 0.2 + (0.6 - 0.2)/2;
shear_centre_pos = 10 + shear_centre_over_c.*chord_dist; %% Needs to be chnaged with actual values

%% Define the locations of each different load where they act
wing_cg_over_c = 0.25 + 0.4894/c_bar_wing;
uc_cg_position = 12.45;
uc_cg_position_over_c = (0.25*chord(uc_cg_position, c0, ctip, b) + 0.95 + 0.4894)./chord(uc_cg_position, c0, ctip, b);
fuel_cg_over_c = 0.2 + (0.6 - 0.2)/2;
engine_location = 4.85 - de_fuselage/2; % m
engine_z_location = 1.2;
engine_cg_location_over_c = 0.189./chord(engine_location, c0, ctip, b);

%% Define the engine thrust load
engine_thrust_load = zeros(1, length(y));
engine_thrust_load(find(y == round(engine_location/ds)*ds)) = Thrust;
% find(y == round(engine_location/ds)*ds)
%% Define the arrrays for torque dist and combined torque
torque = zeros(1, length(y));
comb_torque = zeros(1, length(y));
%% Determine the torque distributions for each different load component for each wing station per unit span
Aero_load_torque = -lift_per_span.*(shear_centre_over_c - x_ac).*chord_dist;
wing_weight_torque = wing_weight_per_span.*(shear_centre_over_c - wing_cg_over_c).*chord_dist;
fuel_weight_torque = fuel_weight_per_span.*(shear_centre_over_c - fuel_cg_over_c).*chord_dist;
engine_weight_torque = engine_weight_load.*(shear_centre_over_c - engine_cg_location_over_c).*chord_dist;
uc_weight_torque = uc_weight_per_span.*(shear_centre_over_c - uc_cg_position_over_c);
engine_thrust_torque = engine_thrust_load.*(-engine_z_location);
%% Determine the torque dists at each span station suing trapezium method
for i=1:length(y)-1
    
    Aero_load_torque(i) = (Aero_load_torque(i) + Aero_load_torque(i+1))*ds/2;

    M_0_torque(i) = (M_0_torque(i) + M_0_torque(i+1))*ds/2;

    wing_weight_torque(i) = (wing_weight_torque(i) + wing_weight_torque(i+1))*ds/2;

    fuel_weight_torque(i) = (fuel_weight_torque(i) + fuel_weight_torque(i+1))*ds/2;

    engine_weight_torque(i) = (engine_weight_torque(i));

    uc_weight_torque(i) = (uc_weight_torque(i) + uc_weight_torque(i+1))*ds/2;

    engine_thrust_torque(i) = (engine_thrust_torque(i));

end


Aero_load_torque(length(y)) = 0;
M_0_torque(length(y)) = 0;
wing_weight_torque(length(y)) = 0;
fuel_weight_torque(length(y)) = 0;
engine_weight_torque(length(y)) = 0;
uc_weight_torque(length(y)) = 0;
engine_thrust_torque(length(y)) = 0;

comb_torque = Aero_load_torque + M_0_torque + wing_weight_torque + fuel_weight_torque + engine_thrust_torque + engine_weight_torque + uc_weight_torque;


%% Find the total torque distribution along the whole span of the wing
for i=1:length(y)
    
    torque(i) = sum(comb_torque(i:end));

end