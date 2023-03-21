function torque = HTTorque(HT, weight, lift, y)
% Takes discretised horizontal tail (HT) and its total weight and lift
% outputs torque (LE down) at each station of HT

index       = length(HT);
ds          = abs(HT(2)-HT(1));

% chordwise positions
x_cg        = 0.4;
x_fspar     = 0.25;
x_rspar     = 0.72;
x_flexAx    = (x_fspar+x_rspar)/2;

% other geometric parameters
fuselage    = 0.854;
HT          = HT + fuselage;
taper       = 0.4;
phi_4       = 10*pi/180;
AR          = 4.5;

% distances for moment arms
tanPhi_FA   = tan(phi_4) - 4/AR * (x_flexAx-0.25)*(1-taper)/(1+taper);
tanPhi_cg   = tan(phi_4) - 4/AR * (x_cg-0.25)*(1-taper)/(1+taper);
tanPhi_4    = 10*pi/180;
flexAx      = HT * tanPhi_FA;
c4          = HT * tanPhi_4;
cg          = HT * tanPhi_cg;

% initialise torque arrays
T_aero = zeros(size(HT));
T_weight = zeros(size(HT));

% find torques by integrating product of sectional forces and moment arms
for i = 1:index
    T_aero(i)   = ds*trapz(lift(i:index) .* (c4(i:index) - flexAx(i)) );
    T_weight(i) = ds*trapz(weight(i:index) .* (flexAx(i) - cg(i:index)) );
end

% combine in single array
torque = T_aero + T_weight;

[~,idx] = min(abs(y-HT));
torque = torque(idx);
end