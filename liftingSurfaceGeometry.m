% Housekeeping
clear
clc

%% Load Polar
load naca65412.mat
load naca64010.mat

%% Main Wing
% Conceptual parameters
wing.span   = 27.672;
wing.C_r    = 4.39;
wing.taper  = 0.4;
wing.sweep  = 0;
wing.polar  = naca65412;
wing.LEfus  = 3.344/2;
wing.TEfus  = 3.344/2;
wing.LEfusC = wing.C_r*(1-(1-wing.taper)*wing.LEfus/(wing.span/2));
wing.TEfusC = wing.C_r*(1-(1-wing.taper)*wing.TEfus/(wing.span/2));

% Spars
wing.fspar          = 0.2;
wing.rspar          = 0.6;

% Find wing geometry parameters
wing.coords = perimeter(wing);
[wing.fsparCoords, wing.rsparCoords]  = spar(wing);

% Wingbox
wing.wingboxCoords  = wingbox(wing);
wing.maxHeight      = 0.12*wing.LEfusC;
wing.minHeight      = 0.12*wing.C_r*wing.taper;
wing.wingboxWidth   = (wing.rspar-wing.fspar);

%% Horizontal Tail
% Conceptual parameters
HT.span     = 9.538;
HT.C_r      = 3.028;
HT.taper    = 0.4;
HT.sweep    = 10;
HT.polar    = naca64010;
HT.LEfus    = 1.1509;
HT.TEfus    = 0.3618;
HT.LEfusC = HT.C_r*(1-(1-HT.taper)*HT.LEfus/(HT.span/2));
HT.TEfusC = HT.C_r*(1-(1-HT.taper)*HT.TEfus/(HT.span/2));

% Spars
HT.fspar  = 0.15;
HT.rspar  = 0.65;

% Find wing geometry parameters
HT.coords = perimeter(HT);
[HT.fsparCoords, HT.rsparCoords]  = spar(HT);

% Wingbox
HT.wingboxCoords  = wingbox(HT);
HT.maxHeight      = 0.12*HT.C_r*(1-(1-HT.taper)*HT.LEfus/(HT.span/2));
HT.minHeight      = 0.12*HT.C_r*HT.taper;
HT.wingboxWidth   = (HT.rspar-HT.fspar);

%% Vertical Tail
% Conceptual parameters
VT.span     = 5.585*2;
VT.C_r      = 4.432;
VT.taper    = 0.4;
VT.sweep    = 35;
VT.polar    = naca64010;
VT.LEfus    = 0;
VT.TEfus    = 0;

% Spars
VT.fspar  = 0.15;
VT.rspar  = 0.58;

% Find wing geometry parameters
VT.coords = perimeter(VT);
[VT.fsparCoords, VT.rsparCoords]  = spar(VT);

%% Wingbox
VT.wingboxCoords  = wingbox(VT);
VT.maxHeight      = 0.12*VT.C_r*(1-(1-VT.taper)*VT.LEfus/VT.span);
VT.minHeight      = 0.12*VT.C_r*VT.taper;

%% Save to .mat
save('liftSurfGeom.mat', 'wing', 'VT', 'HT');