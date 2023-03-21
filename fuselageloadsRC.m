clc
clear all

% The first step is to equilibriate moments about station 0 due to weight
% acting about the CG and reaction forces at the rear and front spar.

% Let the following parameters be reaction forces at the front and rear
% spars respectively, positive up. 

% N.B. This is for inertial loading. Aerodynamic loading is introduced
% with subscripts fa and ra.

R_fi = 0;
R_ri = 0;

%The weight of the plane and centre of gravity has been calculated
%previously

total_mass = 45119; % kg
total_weight= total_mass*9.81; % N


x_cg = 11.5272; % m

% The front spar is at the quarter chord and the rear spar is at 67.5% of
% the chord due to the presence of ailerons. Hence, it is possible to
% determine the location of x_f and x_r, where x_f and x_r are the x
% locations of the front and rear spars where the wing contacts the
% fuselage respectively

x_f = 0;
x_r = 0;

% N.B: x_w is 11.5 m which is the distance from nose to the leading edge of
% the wing at the root

x_w = 11.5; % m

% Henceforth, we need to determine the chord at the root, given that the
% wing 'extends' to the midline

root_chord = 4.0739; % m

% Hence:

x_f = x_w + (root_chord * 0.25); % m 
x_r = x_w + (root_chord * 0.675); % m

% Moments about the nose:
% -(total_weight * x_cg) + (R_f * x_f) + (R_r * x_r) = 0 [1]

% Vertical force equilibrium:
% -total_weight + R_f + R_r = 0 [2]

% Solving [1],[2] simultaneoulsly yields:

R_ri = ((total_weight * x_cg) - (x_f * total_weight))/ (x_r - x_f) % N
R_fi = total_weight - R_ri % N

% Now, assume an input tail aerodynamic load; say 'L', and set to an
% arbitrary value which can later be manipulated

L_kg = 100; % kg
L = L_kg * 9.81; % N

% Asumme that this tail load acts at the quarter-chord of the horizontal
% tail, denoted by x_l. Henceforth:

HT_root_chord = 3.0279; % m
x_HT = 25.0536; % m
x_l = x_HT + (0.25 * HT_root_chord);

% Due to this load at the tail, it is possible to determine R_ra and R_fa
% due to aero load

% N.B. The loading on the tail is defined as downwards
% Moments about the nose:
% (R_f * x_f) + (R_r * x_r) - (L * x_l) = 0 [1]

% Vertical force equilibrium:
% -L + R_f + R_r = 0 [2] 

R_ra = ((L * x_l) - (x_f * L))/ (x_r - x_f) % N
R_fa = L - R_ra % N

% ---------------- 
% Solely consider aero loads:

% Now, we have a simple problem in which three forces on the fuselage of
% the plane generate shear forces and bending moment. The distribution
% across the fuselage can be determined:

% A note on convention: Shear down (L->R) is denoted positive and bending 
% moment which leads to hogging is denoted positive

% At the front spar:

S_f = R_fa; % N
S_r = R_ra; % N
S_L = -L; % N

% From 0 < x < x_f we have 0 shear [S1]
% From x_f < x < x_r we have S_f shear [S2]
% From x_r < x < x_l we have S_f + S_r shear [S3]
% From x_l < x < end we have S_f + S_r + S_L shear, which should be zero
% [S4]

% Define x1,x2,x3,x4 which are discretised points within the ranges above,
% separated by ds

ds=0.1;

S1 = 0; % N
x1=0:ds:x_f;
S2 = S_f; % N
x2=x_f:ds:x_r;
S3 = S_f + S_r; % N
x3=x_r:ds:x_l;
S4 = S_f + S_r + S_L; % N
x4=x_l:ds:27.8599;

% Check S4 is 0. Variable CHECK should return '1'

CHECK = (S4==0)

% Now create a matrix, S, with the shear force from nose to plane end. This
% is coded using poor practice but is fine since the distance of where the
% loads act and number of loads is known. 

S=[0];

for i=1:length(x1)
    S(i)=S1;
end

for j=length(x1):length(x2)+length(x1)
    S(j)=S2;
end

for k=length(x2)+length(x1):length(x1)+length(x2)+length(x3)
    S(k)=S3;
end

for l=length(x1)+length(x2)+length(x3):length(x1)+length(x2)+length(x3)+length(x4)
    S(l)=S4;
end

% Now, work on bending moments using the simple formulation developed by
% RO, MT, MK, RC on 30/01

BM = [0]; 
 
for i = 2:length(S) 
    BM(i) = BM(i-1) + S(i-1)*ds; 
end 

% Now, plot the distributions of shear and bending moment

x=[0];
for i=2:length(S)
    x(i)=x(i-1)+ds;
end

plot(x,S)
title("S v x")
figure
plot(x,BM)
title("BM v x")

% -------------------------------------------------
% Gear Loads:





% Interial Loads:

% First it is necessary to consider interial loads which act on the
% fuselage, the magnitude of each load, the distribution of each load, and 
% the location of each load

% Note: all in kg
%M_Powerplant =	6778.48925;
M_Fuselage =  4368.151897; %Distribute across geoemetry of fuselage
M_HT	=360.0009746; %Distribute across geometry of HT
M_VT	=410.4476968; %Distribute across geometry of VT
M_HandlingGear	=14.10404979; %?
M_FlightControls =	494.9470656; %Point in cockpit
M_FuelSystem = 141.920883; %?
M_Fuel = 14330; %?

M_Fuel_Total = (M_Fuel*(5.22)/18.2) + M_FuelSystem;

M_Avionics	=834.7275343; %?
M_APU=	129.72674; %Distribute APU start to APU end
M_ElectricalSystems	=1037.306605; %?
M_Furnishings	=1819.63202; %Distribute cabin start to cabin end
%M_AntiIcingsystems	=94.02699859;
M_AirConditioning 	=372.759565; %Distribute about location of placement
M_HydraulicSystems	=63.76213039; %?
%M_Aircraftwings	=2666.785041;
%M_MainLandingGear =	605.8813458;
%M_NoseLandingGear	=93.75574989;
M_FrontAirStairs	=100; %Distribute - front door start to front door end
M_BackAirStairs =100; %Distribute - back door start to back door end
M_Water	=400; %Distribute - water tank start to water tank end
M_Instruments	=64.76147073; %Point in cockpit
M_Pilots	=156.9632242; %Distribute - cockpit seat start to cockpit seat 
% end
M_HandBaggage 	=975; %Distribute - cabin start to cabin end
M_Baggage	=1530; %Distribute - hold start to hold end
M_Passengers = (90/93)*7298.789924; %Distribute - cabin start to cabin end
M_Front_Crew = M_Passengers/90;
M_Back_Crew = M_Front_Crew*2;

%N.B. crew and clothing included within passengers and crew luggage
%included within hand baggage

%% Geometry based values, locations in the x direction
Fuselage_length = floor((2785.590043/100)*10)/10;
Cabin_length = round(20.49,1);
Cockpit_length = 4.0;
x_TE_HT = round(26.55590043,1);
HT_root_chord = round(3.02785573,1);
VT_root_chord = round(4.432362376,1);
x_LE_VT = round(21.16403829,1);
x_HandlingGear = round(12.35,1);
x_FlightControls = 1.0;
x_Fuel = 11.5;
x_Stairs = round(12.1173165,1);
x_Water = 24.0;
x_AirConditioning = 3.0;
x_pilot = 2.0;
pilot_seat_length = round(0.91,1);
Cabins_useful_length = round(15.49,1);
x_Cabin_useful_LE = 6.0; 
x_Front_Door_LE = 5.0;
Door_length = 1.0;
Seat_length_standard = round(0.81,1);
x_Back_Seats_LE = x_Cabin_useful_LE+Cabins_useful_length; 
x_Back_Door_LE = x_Back_Seats_LE + Seat_length_standard;



%%  Discretisation of masses, first row is the x location in increments of 0.1m, subsequent rows are detailed below
disc = round([0.0:0.1:Fuselage_length],1);

%Fuselage disc along row 2
disc(2,:) = M_Fuselage/size(disc,2);

%HT disc along row 3
disc(3,(find(disc(1,:) == (x_TE_HT-HT_root_chord)):find(disc(1,:) == (x_TE_HT-HT_root_chord))+(HT_root_chord/0.1))) = M_HT/((HT_root_chord/0.1)+1);

%VT disc along row 4
disc(4,(find(disc(1,:) == (x_LE_VT)):find(disc(1,:) == (x_LE_VT))+(VT_root_chord/0.1))) = M_VT/((VT_root_chord/0.1)+1);

%Handling gear disc along row 5
disc(5,(disc(1,:) == x_HandlingGear)) = M_HandlingGear;

%Flight Controls disc along row 6
disc(6,(disc(1,:) == x_FlightControls)) = M_FlightControls;

%Fuel disc along row 7
disc(7,(disc(1,:) == x_Fuel)) = M_Fuel_Total;

%Avionics disc along row 8
disc(8,(1:(find(disc(1,:) == Cockpit_length)))) = M_Avionics/((4/0.1)+1);

%APU disc along row 9
disc(9,(find(disc(1,:) == x_TE_HT):end)) = M_APU/(size(disc,2)-find(disc(1,:) == x_TE_HT)+1);

%Electrical Systems disc along row 10
disc(10,:) = M_ElectricalSystems/(size(disc,2));

%Furnishings disc along row 11
disc(11,(find(disc(1,:) == x_Cabin_useful_LE):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_Furnishings/((Cabins_useful_length/0.1)+1);

%AirCon disc along row 12
disc(12,(disc(1,:) == x_AirConditioning)) = M_AirConditioning;

%Hydraulic Systems along row 13
disc(13,:) = M_HydraulicSystems/size(disc,2);

%Front Stairs disc alomg row 14
disc(14,(find(disc(1,:) == x_Front_Door_LE):find(disc(1,:) == (x_Front_Door_LE+Door_length)))) = M_FrontAirStairs/((Door_length/0.1)+1);

%Back Stairs disc alomg row 15
disc(15,(find(disc(1,:) == x_Back_Door_LE):find(disc(1,:) == (x_Back_Door_LE+Door_length)))) = M_BackAirStairs/((Door_length/0.1)+1);

%Water disc along row 16
disc(16,(disc(1,:) == x_Water)) = M_Water;

%Instruments disc along row 17
disc(17,(1:(find(disc(1,:) == Cockpit_length)))) = M_Instruments/((4/0.1)+1);

%Pilot disc along row 18
disc(19,(find(disc(1,:) == x_pilot):find(disc(1,:) == pilot_seat_length+x_pilot))) = M_Pilots/((pilot_seat_length/0.1)+1);

%Hand Baggage along row 19
disc(19,(find(disc(1,:) == x_Cabin_useful_LE):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_HandBaggage/((Cabins_useful_length/0.1)+1);

%Baggage along row 20
disc(20,(find(disc(1,:) == x_Cabin_useful_LE):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_Baggage/((Cabins_useful_length/0.1)+1);

%Passengers along row 21
disc(21,(find(disc(1,:) == x_Cabin_useful_LE):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_Passengers/((Cabins_useful_length/0.1)+1);

%Front Crew along row 22
disc(22,(find(disc(1,:) == Cockpit_length):find(disc(1,:) == (Cockpit_length+Seat_length_standard)))) = M_Front_Crew/((Seat_length_standard/0.1)+1);

%Back Crew along row 23
disc(23,(find(disc(1,:) == x_Back_Seats_LE):find(disc(1,:) == (x_Back_Seats_LE+Seat_length_standard)))) = M_Back_Crew/((Seat_length_standard/0.1)+1);

%Total mass and weight
Total_disc_mass = sum(disc((2:end),:),1);
Total_disc_weight = Total_disc_mass*9.81;
Total_weightactingon_fuselage=sum(Total_disc_weight);

Total_moment=sum(Total_disc_weight.*disc(1,:));
Answer_matrix=[Total_moment; Total_weightactingon_fuselage];
Solve_matrix=[x_f x_r; 1 1];
Answer= linsolve(Solve_matrix,Answer_matrix);
R_fi=Answer(1);
R_ri=Answer(2);

%y from discretisation used
y_calc=disc(1,:);

% Reaction forces row 24
Total_disc_weight=Total_disc_weight*-1;
Reactions_i=zeros(1,length(Total_disc_weight));
Reactions_i((y_calc==round(x_f,1)))=R_fi;
Reactions_i((y_calc==round(x_r,1)))=R_ri;

Total_disc_forces=Total_disc_weight+Reactions_i;

%Calculate BM and SF throughout fuselage due to inertial loading
[SF_i,BM_i]=CalculateSFandBM_fuselage(Total_disc_forces,y_calc);

% %Plot due to inertial loading
% figure
% plot(y_calc,SF_i)
% figure
% plot(y_calc,BM_i)

%---------------------------------------
%Combination of aero and inertial loads
% tol = 1;% Tolerance of BM (Nm)
% tol_check = 0;
% L = 0.5;
% 
% while tol_check ~= 1
%     
%     L = L + 0.5;
%         
%     Total_disc_weight2 = Total_disc_mass*-9.81;
%     Total_disc_weight2((y_calc==round(x_l,1)))=Total_disc_weight2((y_calc==round(x_l,1)))+L;
%     Total_weightactingon_fuselage2=sum(Total_disc_weight2);
%     
%     Total_moment2=sum(Total_disc_weight2.*disc(1,:));
%     Answer_matrix=[Total_moment2; Total_weightactingon_fuselage2];
%     Solve_matrix=[x_f x_r; 1 1];
%     Answer= linsolve(Solve_matrix,Answer_matrix);
%     R_f=Answer(1);
%     R_r=Answer(2);
%     
% %     Total_disc_weight2=Total_disc_weight2*-1;
%     Reactions=zeros(1,length(Total_disc_weight));
%     Reactions((y_calc==round(x_f,1)))=R_f;
%     Reactions((y_calc==round(x_r,1)))=R_r;
%     
%     Total_disc_forces2=Total_disc_weight2+Reactions;
%     
%     %Calculate BM and SF throughout fuselage due to inertial loading
%     [SF,M]=CalculateSFandBM_fuselage(Total_disc_forces2,y_calc);
%     
%     
%     if abs(M(end)) < tol 
%         tol_check = 1;
%     end
% 
% end

% %Plot due to inertial loading
% figure
% plot(y_calc,SF)
% figure
% plot(y_calc,M)

% Gear Loads

x_gear_r = round(12.45,1);

tol = 100;% Tolerance of BM (Nm)
tol_check = 0;
L = -0;
n = 2.7;

while tol_check ~= 1
    
    L = L - 10;
        
    Total_disc_weight3 = Total_disc_mass*9.81;
    Total_disc_weight3((y_calc==round(x_gear_r,1)))=Total_disc_weight3((y_calc==round(x_gear_r,1)))-(L+Total_weightactingon_fuselage*n);
    Total_disc_weight3((y_calc==round(x_l,1)))=Total_disc_weight3((y_calc==round(x_l,1)))-L;
    
    
    Total_disc_weight3=Total_disc_weight3*-1;

%     Calculate BM and SF throughout fuselage due to inertial loading
    [SF3,M3]=CalculateSFandBM_fuselage(Total_disc_weight3,y_calc);
    
    
    if abs(M3(end)) < tol 
        tol_check = 1;
    end
%     tol_check = 1;

end

%% Plot due to inertial loading

figure()
hold on
plot(y_calc,SF3)
plot(y_calc,M3)
hold off
grid on


