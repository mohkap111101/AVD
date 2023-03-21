function [SF,BM,disc,L] = FuselageSFandBM(load_case,n,discretization,num)

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
x_l= round(25.8106/num)*num;

% N.B: x_w is 11.5 m which is the distance from nose to the leading edge of
% the wing at the root

x_w = 11.5; % m

% Henceforth, we need to determine the chord at the root, given that the
% wing 'extends' to the midline

root_chord = 4.0739; % m

% Hence:

x_f = x_w + (root_chord * 0.25); % m 
x_r = x_w + (root_chord * 0.675); % m


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
Fuselage_length = floor((2785.590043/100)/num)*num;
Cabin_length = round(20.49/num)*num;
Cockpit_length = round(4.0/num)*num;
x_TE_HT = round(26.55590043/num)*num;
HT_root_chord = round(3.02785573/num)*num;
VT_root_chord = round(4.432362376/num)*num;
x_LE_VT = round(21.16403829/num)*num;
x_HandlingGear = round(12.35/num)*num;
x_FlightControls = round(1.0/num)*num;
x_Fuel = round(11.5/num)*num;
x_Stairs = round(12.1173165/num)*num;
x_Water = round(24.0/num)*num;
x_AirConditioning = round(3.0/num)*num;
x_pilot = round(2.0/num)*num;
pilot_seat_length = round(0.91/num)*num;
Cabins_useful_length = round(15.49/num)*num;
x_Cabin_useful_LE = round(6.0/num)*num;
x_Front_Door_LE = round(5.0/num)*num;
Door_length = round(1.0/num)*num;
Seat_length_standard = round(0.81/num)*num;
x_Back_Seats_LE = x_Cabin_useful_LE+Cabins_useful_length; 
x_Back_Door_LE = x_Back_Seats_LE + Seat_length_standard;



%%  Discretisation of masses, first row is the x location in increments of 'discretization' m, subsequent rows are detailed below
disc = round([0:discretization:Fuselage_length]./num).*num;
%Fuselage disc along row 2
disc(2,:) = M_Fuselage/size(disc,2);

%HT disc along row 3
disc(3,(find(disc(1,:) == (x_TE_HT-HT_root_chord)):find(disc(1,:) == (x_TE_HT-HT_root_chord))+(HT_root_chord/discretization))) = M_HT/((HT_root_chord/discretization)+1);

%VT disc along row 4
disc(4,(find(disc(1,:) == (x_LE_VT)):find(disc(1,:) == (x_LE_VT))+(VT_root_chord/discretization))) = M_VT/((VT_root_chord/discretization)+1);

%Handling gear disc along row 5
disc(5,(disc(1,:) == x_HandlingGear)) = M_HandlingGear;

%Flight Controls disc along row 6
disc(6,(disc(1,:) == x_FlightControls)) = M_FlightControls;

%Fuel disc along row 7
disc(7,(disc(1,:) == x_Fuel)) = M_Fuel_Total;

%Avionics disc along row 8
disc(8,(1:(find(disc(1,:) == Cockpit_length)))) = M_Avionics/((4/discretization)+1);

% %APU disc along row 9
% size(disc)
% find(disc(1,:) == x_TE_HT)
% disc(1,:)
% x_TE_HT
disc(9,(find(disc(1,:) == x_TE_HT):end)) = M_APU/(size(disc,2)-find(disc(1,:) == x_TE_HT)+1);

%Electrical Systems disc along row 10
disc(10,:) = M_ElectricalSystems/(size(disc,2));

%Furnishings disc along row 11
disc(11,(find(disc(1,:) == x_Cabin_useful_LE):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_Furnishings/((Cabins_useful_length/discretization)+1);

%AirCon disc along row 12
disc(12,(disc(1,:) == x_AirConditioning)) = M_AirConditioning;

%Hydraulic Systems along row 13
disc(13,:) = M_HydraulicSystems/size(disc,2);

%Front Stairs disc alomg row 14
disc(14,(find(disc(1,:) == x_Front_Door_LE):find(disc(1,:) == (x_Front_Door_LE+Door_length)))) = M_FrontAirStairs/((Door_length/discretization)+1);

%Back Stairs disc alomg row 15
disc(15,(find(disc(1,:) == x_Back_Door_LE):find(disc(1,:) == (x_Back_Door_LE+Door_length)))) = M_BackAirStairs/((Door_length/discretization)+1);

%Water disc along row 16
disc(16,(disc(1,:) == x_Water)) = M_Water;

%Instruments disc along row 17
disc(17,(1:(find(disc(1,:) == Cockpit_length)))) = M_Instruments/((4/discretization)+1);

%Pilot disc along row 18
disc(19,(find(disc(1,:) == x_pilot):find(disc(1,:) == pilot_seat_length+x_pilot))) = M_Pilots/((pilot_seat_length/discretization)+1);

%Hand Baggage along row 19
disc(19,(find(disc(1,:) == x_Cabin_useful_LE):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_HandBaggage/((Cabins_useful_length/discretization)+1);

%Baggage along row 20
disc(20,(find(disc(1,:) == x_Cabin_useful_LE):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_Baggage/((Cabins_useful_length/discretization)+1);

%Passengers along row 21
disc(21,(find(disc(1,:) == x_Cabin_useful_LE):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_Passengers/((Cabins_useful_length/discretization)+1);

%Front Crew along row 22
disc(22,(find(disc(1,:) == Cockpit_length):find(disc(1,:) == (Cockpit_length+Seat_length_standard)))) = M_Front_Crew/((Seat_length_standard/discretization)+1);

%Back Crew along row 23
disc(23,(find(disc(1,:) == x_Back_Seats_LE):find(disc(1,:) == (x_Back_Seats_LE+Seat_length_standard)))) = M_Back_Crew/((Seat_length_standard/discretization)+1);


%Total mass and weight
Total_disc_mass = sum(disc((2:end),:),1);
Total_disc_weight = Total_disc_mass*9.81;
Total_weightactingon_fuselage=sum(Total_disc_weight);

y_calc=disc(1,:);


if (load_case == 1)

    %Combination of aero and inertial loads
    tol = 1000;% Tolerance of BM (Nm)
    tol_check = 0;
    L = -100;
    
    
    
    while tol_check ~= 1
        
        L = L + 100;
            
        Total_disc_weight = Total_disc_mass*-9.81*n;
        Total_disc_weight(y_calc==x_l)=Total_disc_weight(y_calc==x_l)+L;
        Total_weightactingon_fuselage=sum(Total_disc_weight);
        
        Total_moment=sum(Total_disc_weight.*disc(1,:));
        Answer_matrix=[Total_moment; Total_weightactingon_fuselage];
        Solve_matrix=[x_f x_r; 1 1];
        Answer= linsolve(Solve_matrix,Answer_matrix);
        R_f=Answer(1);
        R_r=Answer(2);
        
%         Total_disc_weight=Total_disc_weight*-1;
        Reactions=zeros(1,length(Total_disc_weight));
        Reactions((y_calc==(round(x_f/num)*num)))=R_f;
        Reactions((y_calc==(round(x_r/num)*num)))=R_r;
        
        
        Total_disc_forces=Total_disc_weight+Reactions;
        
        %Calculate BM and SF throughout fuselage due to inertial loading
        [SF,BM]=CalculateSFandBM_fuselage(Total_disc_forces,y_calc);
        
        
        if abs(BM(end)) < tol 
            tol_check = 1;
        end
        BM(end);

        tol_check = 1;
    end

elseif (load_case == 2)

    x_gear_r = round(12.45/num)*num;
    
    tol = 100;% Tolerance of BM (Nm)
    tol_check = 0;
    L = +5;
    
    while tol_check ~= 1
        
        L = L - 5;
            
        Total_disc_weight = Total_disc_mass*9.81*n;
        Total_disc_weight((y_calc==round(x_gear_r/num)*num))=Total_disc_weight((y_calc==round(x_gear_r/num)*num))-(L+Total_weightactingon_fuselage);
        Total_disc_weight((y_calc==round(x_l/num)*num))=Total_disc_weight((y_calc==round(x_l/num)*num))-L;
        
        
        Total_disc_weight=Total_disc_weight*-1;
    
        %Calculate BM and SF throughout fuselage due to inertial loading
        [SF,BM]=CalculateSFandBM_fuselage(Total_disc_weight,y_calc);
        
        
        if abs(BM(end)) < tol 
            tol_check = 1;
        end
    end

end


