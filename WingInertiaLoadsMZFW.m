function [loads, wing_weight_per_span, fuel_weight_per_span, engine_weight_load, uc_weight_per_span] = WingInertiaLoadsMZFW(y, n, chord_dist)
% Function to determine the distributions of the combined inertial loads 
% along the span of the wing. Combined loading used to calculate the shear
% force and bending moment distributions along the span of the wing.
% Outputs -> y          array of stations along wing span
%            loads      array of combined load at each station
%            SF         array of shear force at each station
%            BM         array of bending moment at each station


load Mohkap_workspace.mat
%% Parameters

% Parameters for loading estimation
s = b/2; % m
f_width = de_fuselage; % m
single_wing_mass = 2667/2; % kg
single_wing_box_span = (b-f_width)/2;
fuel_tank_span = 0.9*single_wing_box_span; % m - should be for only one wing
fuel_capacity = 0; % m^3
engine_location = 4.85 - f_width/2; % m
engine_mass = 2809 + 514; % kg
avg_density = 792; % kgm-^3
uc_length = 1.54; % m
uc_mass = 303; % kg
uc_location_span = 2-f_width/2; % m
uc_mass_in_wing = 65; % kg  - assuming distributed weight along all length

lift_required = MTOW * 9.81 * n; % N
uc_weight_in_wing = uc_mass_in_wing*9.81;
engine_weight = engine_mass * 9.81; % N
fuel_weight_single_wing = fuel_capacity * avg_density * 9.81 % N
single_wing_weight = single_wing_mass * 9.81; % N
total_weight = uc_weight_in_wing + engine_weight + fuel_weight_single_wing + single_wing_weight; % N

%% Discretisation
y = round(y, 1)

%% Determnining Loads
% Define arrays to store the loads for each type of load determined above
% along the span of the wing
wing_weight_per_span = zeros(1, length(y));
fuel_weight_per_span = zeros(1, length(y));
engine_weight_load = zeros(1, length(y));
uc_weight_per_span = zeros(1, length(y));
combined_weight = zeros(1, length((y)));

% Find the indexes of the positions of the different forces along the span
idx_fuel_span_pos = find(y == round(fuel_tank_span, 1));
idx_engine_pos = find(y == round(engine_location, 1));
idx_uc_span_pos = find(y == round(uc_location_span, 1));

% Populate the initialised arrays
wing_weight_per_span(:) = (n*single_wing_weight/Sref) .* chord_dist .* (y(2) - y(1));
fuel_weight_per_span(1:idx_fuel_span_pos) = n*fuel_weight_single_wing/fuel_tank_span .* (y(2) - y(1));
engine_weight_load(idx_engine_pos) = n*engine_weight;
uc_weight_per_span(1:idx_uc_span_pos) = n*uc_weight_in_wing/uc_location_span .* (y(2) - y(1));

combined_weight = (wing_weight_per_span + fuel_weight_per_span + engine_weight_load + uc_weight_per_span);

loads = combined_weight;

