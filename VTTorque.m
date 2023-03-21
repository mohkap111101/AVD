function torque = VTTorque(VT, lift)
% Takes discretised vertical tail (VT) and its total weight and lift
% outputs torque at each station of VT

index       = length(VT);
ds          = abs(VT(2)-VT(1));

% chordwise positions
x_fspar     = 0.25;
x_rspar     = 0.65;
x_flexAx    = (x_fspar+x_rspar)/2;

% other geometric parameters
taper       = 0.4;
phi_4       = 35*pi/180;
AR          = 1.8;

% distances for moment arms
tanPhi_FA   = tan(phi_4) - 4/AR * (x_flexAx-0.25)*(1-taper)/(1+taper);
tanPhi_4    = tan(phi_4);
flexAx      = VT * tanPhi_FA;
c4          = VT * tanPhi_4;

% initialise torque arrays
torque = zeros(size(VT));

% find torques by integrating product of sectional forces and moment arms
for i = 1:index
    torque(i)   = ds*trapz(lift(i:index) .* (c4(i:index) - flexAx(i)) );
end

end