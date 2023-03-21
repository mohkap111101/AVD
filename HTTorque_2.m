function [comb_torque, torque] = HTTorque_2(lift_per_span, HT_weight_per_span, chord_dist, y, x_ac, M_0_torque, c_bar_wing, ctip, c0, b)
% Function to determine the torque distributions and total torque along the
% span of the wing
% Outputs -> y              array of stations along the wing span 
%            comb_torque    array of the torques at each station along the
%                           span
%            torque         array of the total torque at each station
%                           along the span

load Mohkap_workspace.mat

y = round(y, 1);
ds = y(2)-y(1);

%% Define position of shear centre with respect to chord and the fuselage nose

shear_centre_over_c = 0.15 + (0.65 - 0.15)/2;
shear_centre_pos = 10 + shear_centre_over_c.*chord_dist; %% Needs to be chnaged with actual values

%% Define the locations of each different load where they act
HT_cg_over_c = 0.25 + 0.4894/c_bar_tail;

%% Define the arrrays for torque dist and combined torque
torque = zeros(1, length(y));
comb_torque = zeros(1, length(y));
%% Determine the torque distributions for each different load component for each wing station per unit span
Aero_load_torque = ds*lift_per_span.*(shear_centre_over_c - x_ac).*chord_dist;
HT_weight_torque = ds*HT_weight_per_span.*(shear_centre_over_c - HT_cg_over_c).*chord_dist;
%% Determine the torque dists at each span station suing trapezium method
for i=1:length(y)-1
    
    Aero_load_torque(i) = (Aero_load_torque(i) + Aero_load_torque(i+1))*ds/2;

    M_0_torque(i) = (M_0_torque(i)*0.1 + M_0_torque(i+1)*0.1)*ds/2;

    HT_weight_torque(i) = (HT_weight_torque(i) + HT_weight_torque(i+1))*ds/2;

end

Aero_load_torque(end), M_0_torque(end), HT_weight_torque(end) = 0;

comb_torque = Aero_load_torque + M_0_torque + HT_weight_torque;

%% Find the total torque distribution along the whole span of the wing
for i=1:length(y)
    
    torque(i) = sum(comb_torque(i:end));

end