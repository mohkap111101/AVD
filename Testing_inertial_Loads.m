function [SF,BM,disc,L,R_f,R_r] = Testing_inertial_Loads(load_case,n,discretization)

testM = [0 405	176	0	139.3	0	185.6	141	34.5	4.5	9.5	53.5	10];
testL = [0.0	-3973.1	-1726.6	13218.8	-1366.5	-3159.1	-1820.7	-1383.2	-338.4	-44.1	-93.2	784.3	-98.1];
testy = [0	27.9	127	185.4	203.2	294.6	304.8	431.8	508	584.2	660.4	736.6	800.1];
%%


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
x_l= round(25.8106/discretization)*discretization;

% N.B: x_w is 11.5 m which is the distance from nose to the leading edge of
% the wing at the root

x_w = 11.5; % m

% Henceforth, we need to determine the chord at the root, given that the
% wing 'extends' to the midline

root_chord = 4.0739; % m

% Hence:

x_f = x_w + (root_chord * 0.25); % m 
x_r = x_w + (root_chord * 0.675); % m

x_f = round(x_f/discretization)*discretization;
x_r = round(x_r/discretization)*discretization;


%%
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

fuselage_total_mass = M_Fuselage+M_HT+M_VT+M_HandlingGear+M_FlightControls+M_Fuel_Total+M_Avionics+M_APU+M_ElectricalSystems+M_Furnishings+M_AirConditioning+M_HydraulicSystems+M_FrontAirStairs+M_BackAirStairs+M_Water+M_Instruments+M_Pilots+M_HandBaggage+M_Baggage+7298.789924;

%N.B. crew and clothing included within passengers and crew luggage
%included within hand baggage

%% Geometry based values, locations in the x direction
Fuselage_length = floor((2785.590043/100)/discretization)*discretization;
Cabin_length = round(20.49/discretization)*discretization;
Cockpit_length = round(4.0/discretization)*discretization;
x_TE_HT = round(26.55590043/discretization)*discretization;
HT_root_chord = round(3.02785573/discretization)*discretization;
VT_root_chord = round(4.432362376/discretization)*discretization;
x_LE_VT = round(21.16403829/discretization)*discretization;
x_HandlingGear = round(12.35/discretization)*discretization;
x_FlightControls = round(1.0/discretization)*discretization;
x_Fuel = round(11.5/discretization)*discretization;
x_Stairs = round(12.1173165/discretization)*discretization;
x_Water = round(24.0/discretization)*discretization;
x_AirConditioning = round(3.0/discretization)*discretization;
x_pilot = round(2.0/discretization)*discretization;
pilot_seat_length = round(0.91/discretization)*discretization;
Cabins_useful_length = round(15.49/discretization)*discretization;
x_Cabin_useful_LE = round(6.0/discretization)*discretization;
x_Front_Door_LE = round(5.0/discretization)*discretization;
Door_length = round(1.0/discretization)*discretization;
Seat_length_standard = round(0.81/discretization)*discretization;
x_Back_Seats_LE = x_Cabin_useful_LE+Cabins_useful_length;
x_Back_Door_LE = x_Back_Seats_LE + Seat_length_standard;
x_gear_r = round(12.45/discretization)*discretization;


%%  Discretisation of masses, first row is the x location in increments of 0.1m, subsequent rows are detailed below
disc = round([0:discretization:Fuselage_length]./discretization).*discretization;

%Fuselage disc along row 2
% disc(2,(2:end)) = M_Fuselage/(size(disc,2)-1);

disc(2,(2:(find(disc(1,:) == x_f)-1))) = M_Fuselage/(size(disc,2)-3);
disc(2,((find(disc(1,:) == x_f)+1):(find(disc(1,:) == x_r)-1))) = M_Fuselage/(size(disc,2)-3);
disc(2,((find(disc(1,:) == x_r)+1):end)) = M_Fuselage/(size(disc,2)-3);


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
disc(8,(2:(find(disc(1,:) == Cockpit_length)))) = M_Avionics/((Cockpit_length/discretization));

% %APU disc along row 9
disc(9,(find(disc(1,:) == x_TE_HT):end)) = M_APU/(size(disc,2)-find(disc(1,:) == x_TE_HT)+1);

%Electrical Systems disc along row 10
% disc(10,(2:end)) = M_ElectricalSystems/(size(disc,2)-1);


disc(10,(2:(find(disc(1,:) == x_f)-1))) = M_ElectricalSystems/(size(disc,2)-3);
disc(10,((find(disc(1,:) == x_f)+1):(find(disc(1,:) == x_r)-1))) = M_ElectricalSystems/(size(disc,2)-3);
disc(10,((find(disc(1,:) == x_r)+1):end)) = M_ElectricalSystems/(size(disc,2)-3);

%Furnishings disc along row 11
disc(11,(find(disc(1,:) == x_Cabin_useful_LE):(find(disc(1,:) == x_f)-1))) = M_Furnishings/((Cabins_useful_length/discretization)-1);
disc(11,((find(disc(1,:) == x_f)+1):(find(disc(1,:) == x_r)-1))) = M_Furnishings/((Cabins_useful_length/discretization)-1);
disc(11,((find(disc(1,:) == x_r)+1):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_Furnishings/((Cabins_useful_length/discretization)-1);

%AirCon disc along row 12
disc(12,(disc(1,:) == x_AirConditioning)) = M_AirConditioning;

%Hydraulic Systems along row 13
% disc(13,(2:end)) = M_HydraulicSystems/(size(disc,2)-1);


disc(13,(2:(find(disc(1,:) == x_f)-1))) = M_HydraulicSystems/(size(disc,2)-3);
disc(13,((find(disc(1,:) == x_f)+1):(find(disc(1,:) == x_r)-1))) = M_HydraulicSystems/(size(disc,2)-3);
disc(13,((find(disc(1,:) == x_r)+1):end)) = M_HydraulicSystems/(size(disc,2)-3);

%Front Stairs disc alomg row 14
disc(14,(find(disc(1,:) == x_Front_Door_LE):find(disc(1,:) == (x_Front_Door_LE+Door_length)))) = M_FrontAirStairs/((Door_length/discretization)+1);

%Back Stairs disc alomg row 15
disc(15,(find(disc(1,:) == x_Back_Door_LE):find(disc(1,:) == (x_Back_Door_LE+Door_length)))) = M_BackAirStairs/((Door_length/discretization)+1);

%Water disc along row 16
disc(16,(disc(1,:) == x_Water)) = M_Water;

%Instruments disc along row 17
disc(17,(2:(find(disc(1,:) == Cockpit_length)))) = M_Instruments/((Cockpit_length/discretization));

%Pilot disc along row 18
disc(18,(find(disc(1,:) == x_pilot):find(disc(1,:) == pilot_seat_length+x_pilot))) = M_Pilots/((pilot_seat_length/discretization)+1);

%Hand Baggage along row 19
% disc(19,(find(disc(1,:) == x_Cabin_useful_LE):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_HandBaggage/((Cabins_useful_length/discretization)+1);

disc(19,(find(disc(1,:) == x_Cabin_useful_LE):(find(disc(1,:) == x_f)-1))) = M_HandBaggage/((Cabins_useful_length/discretization)-1);
disc(19,((find(disc(1,:) == x_f)+1):(find(disc(1,:) == x_r)-1))) = M_HandBaggage/((Cabins_useful_length/discretization)-1);
disc(19,((find(disc(1,:) == x_r)+1):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_HandBaggage/((Cabins_useful_length/discretization)-1);


%Baggage along row 20
% disc(20,(find(disc(1,:) == x_Cabin_useful_LE):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_Baggage/((Cabins_useful_length/discretization)+1);

disc(20,(find(disc(1,:) == x_Cabin_useful_LE):(find(disc(1,:) == x_f)-1))) = M_Baggage/((Cabins_useful_length/discretization)-1);
disc(20,((find(disc(1,:) == x_f)+1):(find(disc(1,:) == x_r)-1))) = M_Baggage/((Cabins_useful_length/discretization)-1);
disc(20,((find(disc(1,:) == x_r)+1):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_Baggage/((Cabins_useful_length/discretization)-1);


%Passengers along row 21
% disc(21,(find(disc(1,:) == x_Cabin_useful_LE):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_Passengers/((Cabins_useful_length/discretization)+1);

disc(21,(find(disc(1,:) == x_Cabin_useful_LE):(find(disc(1,:) == x_f)-1))) = M_Passengers/((Cabins_useful_length/discretization)-1);
disc(21,((find(disc(1,:) == x_f)+1):(find(disc(1,:) == x_r)-1))) = M_Passengers/((Cabins_useful_length/discretization)-1);
disc(21,((find(disc(1,:) == x_r)+1):find(disc(1,:) == (x_Cabin_useful_LE+Cabins_useful_length)))) = M_Passengers/((Cabins_useful_length/discretization)-1);


%Front Crew along row 22
disc(22,(find(disc(1,:) == Cockpit_length):find(disc(1,:) == (Cockpit_length+Seat_length_standard)))) = M_Front_Crew/((Seat_length_standard/discretization)+1);

%Back Crew along row 23
disc(23,(find(disc(1,:) == x_Back_Seats_LE):find(disc(1,:) == (x_Back_Seats_LE+Seat_length_standard)))) = M_Back_Crew/((Seat_length_standard/discretization)+1);

%%
% %Total mass and weight
% Total_disc_mass = sum(disc((2:end),:),1);
% Total_disc_weight = Total_disc_mass*9.81;
% Total_weightactingon_fuselage=sum(Total_disc_weight);

y_calc=disc(1,:);

%Combination of aero and inertial loads
tol = 100;% Tolerance of BM (Nm)
BM = 9000;
L = 4286.4/9.81;
increment = 15;
count = 0;
itr_limit = 1000;

if (load_case == 1)
    figure()
    hold on
    while abs(BM(end)) > tol

        BM_prev = BM(end);

        %Total mass and weight
        Total_disc_mass_inertial = sum(disc((2:end),:),1);


%         Total_disc_mass_inertial = testM;
%         y_calc = round(testy/discretization)*discretization;
%         x_r = round(294.6/discretization)*discretization;
%         x_f = round(185.4/discretization)*discretization;
%         x_l = round(736.6/discretization)*discretization;
        
        Tail_load_mass = zeros(1,length(Total_disc_mass_inertial));
        Tail_load_mass(y_calc == x_l) = -1*L; 

        Total_disc_mass = Total_disc_mass_inertial + Tail_load_mass;

        Total_load = sum(Total_disc_mass);
        Total_Mass_Moment = y_calc * Total_disc_mass';
%         Total_Mass_Moment = sum(Total_disc_mass.*y_calc);

%         Equilibrium_sum = [Total_Mass_Moment;Total_load];
%         Equilibrium_matrix = [x_f, x_r; 1, 1];
%         Reactions = Equilibrium_matrix\Equilibrium_sum;
%         R_f = Reactions(1);
%         R_r = Reactions(2);
        Answer_matrix=[Total_Mass_Moment; Total_load];
        Solve_matrix=[x_f x_r; 1 1];
        Answer= linsolve(Solve_matrix,Answer_matrix);
        R_f=Answer(1);
        R_r=Answer(2);

        Inertial_load_disc = Total_disc_mass_inertial*-9.81*n;

        Aero_load_disc = zeros(1,length(Inertial_load_disc));
        Aero_load_disc(y_calc == x_l) = L*9.81*n;
        Aero_load_disc(y_calc == x_f) = R_f*9.81*n;
        Aero_load_disc(y_calc == x_r) = R_r*9.81*n;

        Total_disc_load = Inertial_load_disc + Aero_load_disc;
    
        %Calculate BM and SF throughout fuselage due to inertial loading
        [SF,BM]=SF_and_BM_daquing(Total_disc_load,y_calc);
       
        plot(y_calc,BM)
        plot(y_calc,SF)
        pause(0.01)
        BM_current = BM(end);

        if abs(BM_current) > abs(BM_prev)
            increment = -1 * increment;
        end

        L = L + increment;

        count = count + 1;
        if (count > itr_limit)
            disp(["Did not converge in,", itr_limit, "iterations"])
            break
        end
    
    end
    hold off
    disp("CONVERGED!!!")
elseif (load_case == 2)
    figure()
    hold on
    while abs(BM(end)) > tol

        BM_prev = BM(end);

        %Total mass and weight
        Total_disc_mass_inertial = sum(disc((2:end),:),1);

        Tail_load_mass = zeros(1,length(Total_disc_mass_inertial));
        Tail_load_mass(y_calc == x_l) = L; 

        Inertial_load_disc = Total_disc_mass_inertial.*-9.81.*n;

        GearLoad = zeros(1,length(Inertial_load_disc));
        GearLoad(y_calc==x_gear_r) = (-n*sum(Total_disc_mass_inertial)+L)*-9.81;
        Tail_Load = Tail_load_mass.*9.81;

        Total_load_disc = Tail_Load + GearLoad + Inertial_load_disc;

%%% ---------------

%         Total_disc_weight = Total_disc_mass.*-9.81*n;
% 
%         Total_Load=abs(sum(Total_disc_weight));
%         
%         Total_disc_weight(y_calc==x_l)=Total_disc_weight(y_calc==x_l) + L;
%         Total_disc_weight(y_calc==x_gear_r) = Total_disc_weight(y_calc==x_gear_r) - (L - Total_Load);
% 
%         Actual_gear_load_factor = 1 + (((-1*(Total_Load + L))+Total_Load)/(-1*Total_Load));
%                
%         Total_Moment = y_calc * Total_disc_weight';
        
    
        %Calculate BM and SF throughout fuselage due to inertial loading
        [SF,BM]=CalculateSFandBM_fuselage(Total_load_disc,y_calc);
        plot(y_calc,BM)
% %         plot(y_calc,SF)
        pause(0.001)
        BM_current = BM(end);

        if abs(BM_current) > abs(BM_prev)
            increment = -1 * increment;
        end

        L = L - increment;

        count = count + 1;
        if (count > itr_limit)
            disp(["Did not converge in,", itr_limit, "iterations"])
            break
        end
    
    end
    hold off
    disp("CONVERGED!!!")
end

L = L*9.81;
R_f = 0;
R_r = 0;

%%
% [testSF,testBM,testBM2]=SF_and_BM_daquing(testL,testy);
% [test2SF,test2BM]=CalculateSFandBM_fuselage(testL,testy);
% % [test3SF,test3BM,test3BM2]=SF_and_BM_daquing(Total_Load_disc,y_calc);
% 
% 
% figure()
% hold on
% plot(testy,testBM,'r');
% plot(testy,test2BM,'b');
% plot(testy,testSF,'g')
% hold off
% grid on


figure()
hold on
plot(y_calc,BM,'r');
hold off
grid on


