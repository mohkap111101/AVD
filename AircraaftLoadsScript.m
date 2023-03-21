% AVD Preliminary Airframe Design Script
% Characterisation of loads


clear;
clc;
close all;

%% Load the variables from workspace file
load Mohkap_workspace.mat

%% V-n diagram for the aircraft

max_lim_load = 2.5; %these are values of load factor
min_lim_load = -1;

max_ult_load = 1.5*max_lim_load;
min_ult_load = 1.5*min_lim_load;

Vc = V_c2; %TAS
Va = sqrt((max_lim_load*MTOW*9.81)/(0.5 * rho_leg2*Sref*Cl_max_aircraft)); %TAS
VH = sqrt((min_lim_load*MTOW*9.81)/(0.5 * rho_leg2*Sref*Cl_min_aircraft)); %TAS

VD = (1/0.8)*Vc; %TAS
V_s = 111.9; %TAS

Va_ult = sqrt((max_ult_load*MTOW*9.81)/(0.5 * rho_leg2*Sref*Cl_max_aircraft)); %TAS
VH_ult = sqrt((min_ult_load*MTOW*9.81)/(0.5 * rho_leg2*Sref*Cl_min_aircraft)); %TAS

% Values of angle of attack and Cl at Va different velocities

Cl_Va = Cl_max_aircraft;
Cl_VD = (max_lim_load * MTOW * 9.81)/(0.5 * rho_leg2 * VD^2 * Sref);

alpha_Va = Cl_Va/a_wing + alpha_0_wing*pi/180;
alpha_VD = Cl_VD/a_wing + alpha_0_wing*pi/80;

% Finding the pithcing moments generated

Cm_0 = -0.075;
Cm_alpha = -0.14132;

M0_Va = (Cm_alpha * alpha_Va + Cm_0) * 0.5 * rho_leg2*Va^2 * Sref * c_bar_wing;
M0_VD = (Cm_alpha * alpha_VD + Cm_0) * 0.5 * rho_leg2 * VD^2 *Sref * c_bar_wing;

% Typipcal values for gust velocities and delta_n 

mu = (2*MTOW*9.81/Sref)/(rho0*9.81*c_bar_wing*a_total);
K = (0.88*mu)/(5.3 + mu);

% Finding VB
U_ref = interp1([15000 60000], [13.41 6.36], h2);

VB = V_s*sqrt(rho_leg2/rho0)*(1 + (K*U_ref*3.281*Vc*1.944*a_total)/(498*MTOW*9.81/Sref * 0.02))^(1/2); % EAS

Ude_VB = 20; % EAS 
Ude_Vc = 15.2; % EAS
Ude_VD = 7.6; %EAS

U_VB = K*Ude_VB * sqrt(rho0/rho0); %TAS
U_Vc = K*Ude_Vc  * sqrt(rho0/rho0); %TAS
U_VD = K*Ude_VD  * sqrt(rho0/rho0); %TAS

delta_n_VB = (rho_leg2*U_VB*VB* sqrt(rho0/rho_leg2) * a_total)/(2*MTOW*9.81/Sref);
delta_n_Vc = (rho_leg2*U_Vc*Vc*a_total)/(2*MTOW*9.81/Sref);
delta_n_VD = (rho_leg2*U_VD*VD*a_total)/(2*MTOW*9.81/Sref);

V_from_0_to_Va = linspace(0, Va, 1000);
n_from_0_to_max = 0.5*rho_leg2*V_from_0_to_Va.^2*Sref*Cl_max_aircraft./(MTOW*9.81);

V_from_0_to_VH = linspace(0, VH, 1000);
n_from_0_to_min = 0.5*rho_leg2*V_from_0_to_VH.^2*Sref*Cl_min_aircraft./(MTOW*9.81);

V_from_Va_to_Va_ult = linspace(Va, Va_ult, 1000);
n_from_max_to_max_ult = 0.5*rho_leg2*V_from_Va_to_Va_ult.^2*Sref*Cl_max_aircraft./(MTOW*9.81);

V_from_VH_to_VH_ult = linspace(VH, VH_ult, 1000);
n_from_min_to_min_ult = 0.5*rho_leg2*V_from_VH_to_VH_ult.^2*Sref*Cl_min_aircraft./(MTOW*9.81);

figure()

% Plotting the manouvre diagram for no gust flight up to limit load
p1 = plot(V_from_0_to_Va.*sqrt(rho_leg2/rho0), n_from_0_to_max, "k", "LineWidth", 2, "DisplayName", "Limit load envelope");
hold on
plot(V_from_0_to_VH.*sqrt(rho_leg2/rho0), n_from_0_to_min, "k", "LineWidth", 2)
plot([Va*sqrt(rho_leg2/rho0) VD*sqrt(rho_leg2/rho0)], [max_lim_load max_lim_load], "k", "LineWidth", 2);
plot([VH*sqrt(rho_leg2/rho0) Vc*sqrt(rho_leg2/rho0)], [min_lim_load min_lim_load], "k", "LineWidth", 2);
plot([VD*sqrt(rho_leg2/rho0) VD*sqrt(rho_leg2/rho0)], [max_lim_load 0], "k", "LineWidth", 2)
plot([Vc*sqrt(rho_leg2/rho0) VD*sqrt(rho_leg2/rho0)], [min_lim_load 0], "k", "LineWidth", 2)
plot([0 0], [-1000 1000], "k", "LineWidth", 1.25);
plot([-1000 1000], [0 0], "k", "LineWidth", 1.25);

% Plotting the manouvre diagram for no gust flight up to limit load
p2 = plot(V_from_Va_to_Va_ult.*sqrt(rho_leg2/rho0), n_from_max_to_max_ult, "r", "LineWidth", 2, "DisplayName", "Ultimate load envelope");
plot(V_from_VH_to_VH_ult.*sqrt(rho_leg2/rho0), n_from_min_to_min_ult, "r", "LineWidth", 2)
plot([Va_ult*sqrt(rho_leg2/rho0), VD*sqrt(rho_leg2/rho0)], [max_ult_load, max_ult_load], "r", "LineWidth", 2)
plot([VH_ult*sqrt(rho_leg2/rho0) Vc*sqrt(rho_leg2/rho0)], [min_ult_load min_ult_load], "r", "LineWidth", 2)
plot([VD*sqrt(rho_leg2/rho0) VD*sqrt(rho_leg2/rho0)], [max_ult_load max_lim_load], "r", "LineWidth", 2)
plot([Vc*sqrt(rho_leg2/rho0) VD*sqrt(rho_leg2/rho0)], [min_ult_load 0], "r", "LineWidth", 2)

%Plotting the gust envelope
%plot([0, VB], [1 1+delta_n_VB], "blue", "LineWidth", 2);
p3 = plot([VB, Vc*sqrt(rho_leg2/rho0)], [1+delta_n_VB 1+delta_n_Vc], "blue", "LineWidth", 2, "DisplayName", "Gust envelope");
plot([VB, Vc*sqrt(rho_leg2/rho0)], [1-delta_n_VB 1-delta_n_Vc], "blue", "LineWidth", 2)
plot([Vc*sqrt(rho_leg2/rho0), VD*sqrt(rho_leg2/rho0)], [1+delta_n_Vc 1+delta_n_VD], "blue", "LineWidth", 2)
plot([Vc*sqrt(rho_leg2/rho0), VD*sqrt(rho_leg2/rho0)], [1-delta_n_Vc 1-delta_n_VD], "blue", "LineWidth", 2)
plot([-1000 1000], [1, 1], "--k", "LineWidth", 1.5);
plot([VD*sqrt(rho_leg2/rho0), VD*sqrt(rho_leg2/rho0)], [1+delta_n_VD 1-delta_n_VD], "blue", "LineWidth", 2)

%Plot velocity indication lines
plot([Va*sqrt(rho_leg2/rho0) Va*sqrt(rho_leg2/rho0)], [-1000, 1000], "-.k", "LineWidth", 1.2);
text(Va*sqrt(rho_leg2/rho0)-5, 3.85, "VA", "Interpreter","latex", "FontSize", 10)
plot([VB VB], [-1000, 1000], "-.k", "LineWidth", 1.2);
text(VB+2, 3.85, "VB", "Interpreter","latex", "FontSize", 10)
plot([Vc*sqrt(rho_leg2/rho0) Vc*sqrt(rho_leg2/rho0)], [-1000, 1000], "-.k", "LineWidth", 1.2);
text(Vc*sqrt(rho_leg2/rho0)-5, 3.85, "VC", "Interpreter","latex", "FontSize", 10)

plot([VD*sqrt(rho_leg2/rho0) VD*sqrt(rho_leg2/rho0)], [-1000, 1000], "-.k", "LineWidth", 1.2);
text(VD*sqrt(rho_leg2/rho0)-5, 3.85, "VD", "Interpreter","latex", "FontSize", 10)
plot([V_s*sqrt(rho_leg2/rho0) V_s*sqrt(rho_leg2/rho0)], [-1000, 1000], "-.k", "LineWidth", 1.2);
text(V_s*sqrt(rho_leg2/rho0)-5, 3.85, "$V_s$", "Interpreter","latex", "FontSize", 10)



%Plot Individual gust points using a line
p4 = plot([0, VB], [1 1+delta_n_VB], "--blue", "LineWidth", 1, "DisplayName", "Gust at VB");
p5 = plot([0, Vc*sqrt(rho_leg2/rho0)], [1 1+delta_n_Vc], "-.blue", "LineWidth", 1, "DisplayName", "Gust at VC");
p6 = plot([0, VD*sqrt(rho_leg2/rho0)], [1 1+delta_n_VD], ":blue", "LineWidth", 1, "DisplayName", "Gust at VD");
plot([0, VB], [1 1-delta_n_VB], "--blue", "LineWidth", 1)
plot([0, Vc*sqrt(rho_leg2/rho0)], [1 1-delta_n_Vc], "-.blue", "LineWidth", 1)
plot([0, VD*sqrt(rho_leg2/rho0)], [1 1-delta_n_VD], ":blue", "LineWidth", 1)

%*sqrt(rho_leg2/rho0)

legend([p1, p2, p3, p4 , p5, p6], "Interpreter","latex", "FontSize", 16)

grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

xlabel("Airspeed (EAS), m/s", "Interpreter","latex", "FontSize", 16)
ylabel("Load factor, n", "Interpreter","latex", "FontSize", 16)

xlim([0, 200])
ylim([-2 4])
hold off

%% Loads on the aircraft structures - Inertial and Aero
% Discretisations
ds = 0.1;
magpie.wing.disc = 0:ds:(b-de_fuselage)/2;

% Wing inertial loads
% Wing aero loads
% Wing combined loads 
% Fuselage inertial loads
% Fuselage aero loads
% Fuselage combined loads
% Tail aero loads
% Tail inertial loads
% Tail combined loads

% Define stations - rib spacing - ds

%% Wing inertial loads
% weight

chord_dist = chord(magpie.wing.disc + de_fuselage/2, c0, c_tip, b);

n = 1;
[magpie.wing.inertial.load.n_1, wing_weight_per_span, fuel_weight_per_span, engine_weight_load, uc_weight_per_span] = WingInertiaLoads(magpie.wing.disc, n, chord_dist);
[magpie.wing.inertial.SF.n_1, magpie.wing.inertial.BM.n_1] = CalculateSFandBM(magpie.wing.inertial.load.n_1, magpie.wing.disc);

% landing
[magpie.wing.gear.SF, magpie.wing.gear.BM] = wingLanding(magpie.wing.disc, MTOW);
magpie.wing.landing.SF = 2.7 * magpie.wing.inertial.SF.n_1 - magpie.wing.gear.SF;
magpie.wing.landing.BM = 2.7 * magpie.wing.inertial.BM.n_1 - magpie.wing.gear.BM;

figure()
plot(magpie.wing.disc, -magpie.wing.inertial.load.n_1, "-k", "LineWidth", 2)
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Inertial Load, N", "Interpreter","latex", "FontSize", 16)
grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")
xlim([0 b/2])

figure()
plot(magpie.wing.disc, magpie.wing.landing.SF, "-ro", "MarkerFaceColor", "r", "MarkerSize", 3)
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Landing Shear Force, N", "Interpreter","latex", "FontSize", 16)
xlim([0 b/2])

figure()
plot(magpie.wing.disc, magpie.wing.landing.BM, "-ro", "MarkerFaceColor", "r", "MarkerSize", 3)
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Landing Bending Moment, Nm", "Interpreter","latex", "FontSize", 16)
xlim([0 b/2])

grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure()
plot(magpie.wing.disc, magpie.wing.inertial.SF.n_1, "-ro", "MarkerFaceColor", "r", "MarkerSize", 3);
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Inertial Shear force, N", "Interpreter","latex", "FontSize", 16)

xlim([0 b/2])

grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure()
plot(magpie.wing.disc, magpie.wing.inertial.BM.n_1, "-ro", "MarkerFaceColor", "r", "MarkerSize", 3)
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Inertial Bending Moment, Nm", "Interpreter","latex", "FontSize", 16)

xlim([0 b/2])

grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

% Inertial loading with load factor 3.75, and at MZFW

[magpie.wing.inertial.load.n_max_ult, wing_weight_per_span_n_ult, fuel_weight_per_span_n_ult, engine_weight_load_n_ult, uc_weight_per_span_n_ult] = WingInertiaLoads(magpie.wing.disc,max_ult_load, chord_dist);
[magpie.wing.inertial.SF.n_max_ult, magpie.wing.inertial.BM.n_max_ult] = CalculateSFandBM(magpie.wing.inertial.load.n_max_ult, magpie.wing.disc);


[magpie.wing.inertialMZFW.load.n_max_ult, wing_weight_per_span_n_ult, fuel_weight_per_span_n_ult, engine_weight_load_n_ult, uc_weight_per_span_n_ult] = WingInertiaLoadsMZFW(magpie.wing.disc,max_ult_load, chord_dist);
[magpie.wing.inertialMZFW.SF.n_max_ult, magpie.wing.inertialMZFW.BM.n_max_ult] = CalculateSFandBM(magpie.wing.inertialMZFW.load.n_max_ult, magpie.wing.disc);

figure()
plot(magpie.wing.disc, magpie.wing.inertialMZFW.SF.n_max_ult)
hold on
plot(magpie.wing.disc, magpie.wing.inertial.SF.n_max_ult)

%% Wing aero loads

% lift distribution
magpie.wing.aero.load.n_1   = wingAero(magpie.wing.disc, MTOW);
magpie.wing.aero.load.n_max = max_ult_load * magpie.wing.aero.load.n_1;
magpie.wing.aero.load.n_min = min_ult_load * magpie.wing.aero.load.n_1;

% zero-lift pitching moment distribution
magpie.wing.pitchingMom.V_a = wingMom(magpie.wing.disc, M0_Va);
magpie.wing.pitchingMom.V_d = wingMom(magpie.wing.disc, M0_VD);

[magpie.wing.aero.SF.n_1, magpie.wing.aero.BM.n_1] = CalculateSFandBM(magpie.wing.aero.load.n_1, magpie.wing.disc);

figure()
plot(magpie.wing.disc, magpie.wing.aero.load.n_1, "-k", "LineWidth", 2)
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Aero load, N", "Interpreter","latex", "FontSize", 16)

xlim([0 b/2])

grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure()
plot(magpie.wing.disc, magpie.wing.aero.SF.n_1, "-bo", "MarkerFaceColor", "r", "MarkerSize", 3);
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Aero Shear force, N", "Interpreter","latex", "FontSize", 16)

xlim([0 b/2])

grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure()
plot(magpie.wing.disc, magpie.wing.aero.BM.n_1, "-bo", "MarkerFaceColor", "r", "MarkerSize", 3)
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Aero Bending Moment, Nm", "Interpreter","latex", "FontSize", 16)

xlim([0 b/2])

grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

%% Wing Combined SF and BM

%, wing_weight_per_span, fuel_weight_per_span, engine_weight_load, uc_weight_per_span] = WingInertiaLoads(magpie.wing.disc, n, chord_dist);
[magpie.wing.inertial.load.n_max_ult, wing_weight_per_span_n_ult, fuel_weight_per_span_n_ult, engine_weight_load_n_ult, uc_weight_per_span_n_ult] = WingInertiaLoads(magpie.wing.disc,max_ult_load, chord_dist);

magpie.wing.combined.load.n_1 = magpie.wing.aero.load.n_1 - magpie.wing.inertial.load.n_1;
[magpie.wing.combined.SF.n_1, magpie.wing.combined.BM.n_1] = CalculateSFandBM(magpie.wing.combined.load.n_1, magpie.wing.disc);

%magpie.wing.combined.BM.n_1 = magpie.wing.aero.BM.n_1 - magpie.wing.inertial.BM.n_1;
%magpie.wing.combined.SF.n_1 = magpie.wing.aero.SF.n_1 - magpie.wing.inertial.SF.n_1;


magpie.wing.combined.load.n_max = magpie.wing.aero.load.n_max - magpie.wing.inertial.load.n_max_ult;
magpie.wing.combinedMZFW.load.n_max = magpie.wing.aero.load.n_max - magpie.wing.inertialMZFW.load.n_max_ult;

[magpie.wing.inertialMZFW.load.n_1, ~, ~, ~, ~] = WingInertiaLoadsMZFW(magpie.wing.disc,1, chord_dist);
[magpie.wing.inertialMZFW.SF.n_1, magpie.wing.inertialMZFW.BM.n_1] = CalculateSFandBM(magpie.wing.inertialMZFW.load.n_1, magpie.wing.disc);

[magpie.wing.combined.SF.n_max, magpie.wing.combined.BM.n_max] = CalculateSFandBM(magpie.wing.combined.load.n_max, magpie.wing.disc);
[magpie.wing.combinedMZFW.SF.n_max, magpie.wing.combinedMZFW.BM.n_max] = CalculateSFandBM(magpie.wing.combinedMZFW.load.n_max, magpie.wing.disc);

figure()
plot(magpie.wing.disc, magpie.wing.combinedMZFW.SF.n_max)
hold on
plot(magpie.wing.disc, magpie.wing.combined.SF.n_max)
hold off

figure()

plot(magpie.wing.disc, magpie.wing.combined.BM.n_max)
hold off

figure()
plot(magpie.wing.disc, magpie.wing.combined.SF.n_1)
hold on
plot(magpie.wing.disc, magpie.wing.aero.SF.n_1)
hold off

figure()
plot(magpie.wing.disc, magpie.wing.combined.BM.n_1)
hold on
plot(magpie.wing.disc, magpie.wing.aero.BM.n_1)
hold off

figure()
plot(magpie.wing.disc, magpie.wing.combined.SF.n_1, "-bo", "MarkerFaceColor", "r", "MarkerSize", 3);
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Combined Shear force, N", "Interpreter","latex", "FontSize", 16)

xlim([0 b/2])

grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure()
plot(magpie.wing.disc, magpie.wing.combined.BM.n_1, "-bo", "MarkerFaceColor", "r", "MarkerSize", 3)
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Combined Bending Moment, Nm", "Interpreter","latex", "FontSize", 16)

xlim([0 b/2])

grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure()
plot(magpie.wing.disc, magpie.wing.combined.load.n_max, "-bo", "MarkerFaceColor", "r", "MarkerSize", 3)
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Combined ultimate load, Nm", "Interpreter","latex", "FontSize", 16)

xlim([0 b/2])

grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")


[magpie.wing.gearMZFW.SF, magpie.wing.gearMZFW.BM] = wingLanding(magpie.wing.disc, MTOW - 7.1478e+03);
magpie.wing.landingMZFW.SF = 2.7 * magpie.wing.inertialMZFW.SF.n_1 - magpie.wing.gearMZFW.SF;
magpie.wing.landingMZFW.BM = 2.7 * magpie.wing.inertialMZFW.BM.n_1 - magpie.wing.gearMZFW.BM;

%% Plots for the report
figure()
p1 = plot(magpie.wing.disc, magpie.wing.combined.SF.n_max, "-k", "LineWidth", 2, "DisplayName","Loading at $V_A$ and $V_D$");
hold on
p2 = plot(magpie.wing.disc, magpie.wing.combinedMZFW.SF.n_max, "-.k", "LineWidth", 2, "DisplayName","Loading at $V_A$ and $V_D$ ZFW");
p3 = plot(magpie.wing.disc, magpie.wing.landing.SF, "-r", "LineWidth", 2, "DisplayName","Landing nose off");
p4 = plot(magpie.wing.disc, magpie.wing.landingMZFW.SF, "--r", "LineWidth", 2, "DisplayName","Landing nose off ZFW");
plot([-1000 1000], [0 0], "-k", "LineWidth", 0.5)
hold off
grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")
legend([p1 p2 p3 p4], "Interpreter","latex", "FontSize", 16)
xlabel("Station along wing (m)", "Interpreter","latex", "FontSize", 16)
ylabel("Shear Force (N)", "Interpreter","latex", "FontSize", 16)
xlim([0 13])


figure()
p1 = plot(magpie.wing.disc, magpie.wing.combined.BM.n_max, "-k", "LineWidth", 2, "DisplayName","Loading at $V_A$ and $V_D$");
hold on
p2 = plot(magpie.wing.disc, magpie.wing.combinedMZFW.BM.n_max, "-.k", "LineWidth", 2, "DisplayName","Loading at $V_A$ and $V_D$ ZFW");
p3 = plot(magpie.wing.disc, magpie.wing.landing.BM, "-r", "LineWidth", 2, "DisplayName","Landing nose off");
p4 = plot(magpie.wing.disc, magpie.wing.landingMZFW.BM, "--r", "LineWidth", 2, "DisplayName","Landing nose off ZFW");
plot([-1000 1000], [0 0], "-k", "LineWidth", 0.5)
hold off
grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")
legend([p1 p2 p3 p4], "Interpreter","latex", "FontSize", 16)
xlabel("Station along wing (m)", "Interpreter","latex", "FontSize", 16)
ylabel("Bending Moment (Nm)", "Interpreter","latex", "FontSize", 16)
xlim([0 13])

%% Wing torque 

V = Va;
chord_dist = chord(magpie.wing.disc + de_fuselage/2, c0, c_tip, b);
Thrust = 143200/2;
T = 1;

[comb_torque_n_ult_vd, torque_n_ult_vd] = WingTorque(magpie.wing.aero.load.n_max, wing_weight_per_span_n_ult, fuel_weight_per_span_n_ult, engine_weight_load_n_ult, uc_weight_per_span_n_ult, chord_dist, magpie.wing.disc, 0.25, magpie.wing.pitchingMom.V_d, Thrust, c_bar_wing, c_tip, c0, b);
[comb_torque_n_ult_va, torque_n_ult_va] = WingTorque(magpie.wing.aero.load.n_max, wing_weight_per_span_n_ult, fuel_weight_per_span_n_ult, engine_weight_load_n_ult, uc_weight_per_span_n_ult, chord_dist, magpie.wing.disc, 0.25, magpie.wing.pitchingMom.V_a, Thrust, c_bar_wing, c_tip, c0, b);
[comb_torque, torque] = WingTorque(magpie.wing.aero.load.n_1, wing_weight_per_span, fuel_weight_per_span, engine_weight_load, uc_weight_per_span, chord_dist, magpie.wing.disc, 0.25, magpie.wing.pitchingMom.V_d, Thrust, c_bar_wing, c_tip, c0, b);
[comb_torque_landing, torque_landing] = WingTorqueLanding(wing_weight_per_span_n_ult, fuel_weight_per_span_n_ult, engine_weight_load_n_ult, chord_dist, magpie.wing.disc, 0.25, 0.3*Thrust, c_bar_wing, c_tip, c0, b);

final_torque = zeros(1,length(magpie.wing.disc));

for i = 1:length(magpie.wing.disc)
    max_torque = 0;
    if(abs(torque_landing(i)) >= abs(torque_n_ult_vd(i)))
        final_torque(i) = torque_landing(i);
    else
        final_torque(i) = torque_n_ult_vd(i);
    end
end

figure()
p1 = plot(magpie.wing.disc, torque_n_ult_va, "-k", "LineWidth", 2, "DisplayName","Loading at $V_A$");

hold on

p2 = plot(magpie.wing.disc, torque_n_ult_va, "-.k", "LineWidth", 2, "DisplayName","Loading at $V_A$ ZFW");
p3 = plot(magpie.wing.disc, torque_n_ult_vd, "-r", "LineWidth", 2, "DisplayName","Loading at $V_D$");
p4 = plot(magpie.wing.disc, torque_n_ult_vd, "-.r", "LineWidth", 2, "DisplayName","Loading at $V_D$ ZFW");
p5 = plot(magpie.wing.disc, torque_landing, "-b", "LineWidth", 2, "DisplayName","Landing nose off");
p6 = plot(magpie.wing.disc, torque_landing, "-.b", "LineWidth", 2, "DisplayName","Landing nose off ZFW");
p7 = plot(magpie.wing.disc, final_torque, "-.g", "LineWidth", 2, "DisplayName","Constraining torque");

plot([-1000 1000], [0 0], "-k", "LineWidth", 0.5)
grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")
legend([p1 p2 p3 p4 p5 p6 p7], "Interpreter","latex", "FontSize", 16)
xlabel("Station along wing (m)", "Interpreter","latex", "FontSize", 16)
ylabel("Torque (Nm)", "Interpreter","latex", "FontSize", 16)
xlim([0 13])


%% Tail loads
% tail weight

load liftSurfGeom.mat
t = 1;

weight_HT_n_max = 360/2*3.75; % CHANGE
weight_HT_n_min = 360/2*-1; % CHANGE

ds = 0.01;

% discretisations
magpie.HT.disc = 0: ds : b_tail/2 - HT.TEfus; % modify numbers
magpie.VT.disc = 0: ds : b_vtail/2; % modify numbers

% key parameters
lift_HT_n_max   = 4.3e3 * 3.75; % CHANGE
lift_HT_n_min   = -4000/2; % CHANGE
lift_VT         = 28685;
lift_HT_landing = 33000;

% horizontal tail max load factor
magpie.HT.lift_landing.n_max        = HTAero(magpie.HT.disc, lift_HT_landing);
magpie.HT.lift.n_max        = HTAero(magpie.HT.disc, lift_HT_n_max);
% magpie.HT.weight.n_max      = HTweight(magpie.HT.disc, weight_HT_n_max);

weight_total = 360*9.81/2;
ht  = linspace(0,HT.span/2, 1000);
w_r = 2 * weight_total/((1+HT.taper)*HT.span);
weight  = w_r * (1 - (1-HT.taper) * (HT.TEfus + ht)/(HT.span/2));
lift_landing = HTAero(ht, lift_HT_landing);
lift2 = (HTAero(ht,lift_HT_n_max));
magpie.HT.torque.n_max      = HTTorque1(magpie.HT.disc, ...
                                weight, lift2);


magpie.HT.torqueLanding.n_max      = HTTorqueLanding(magpie.HT.disc, ...
                                weight, lift_landing);


%magpie.HT.combined.n_max    = magpie.HT.lift.n_max - magpie.HT.weight.n_max; 

%[magpie.HT.SF.n_max, magpie.HT.BM.n_max]    = CalculateSFandBM(...
 %                                               magpie.HT.combined.n_max,...
  %                                              magpie.HT.disc);


%horizontal tail min load factor
% magpie.HT.lift_landing.n_min        = HTAero(magpie.HT.disc, lift_HT_n_min);
% magpie.HT.weight.n_min      = HTweight(HT, weight_HT_n_min);
% magpie.HT.torque.n_min      = HTTorque(magpie.HT.disc, ...
%                                 magpie.HT.weight.n_min, magpie.HT.lift.n_min);
%magpie.HT.combined.n_min    = magpie.HT.lift.n_min - magpie.HT.weight.n_min;

%[magpie.HT.SF.n_min, magpie.HT.BM.n_min]    = CalculateSFandBM(...
 %                                               magpie.HT.combined.n_min,...
  %                                              magpie.HT.disc);


% vertical tail OEI
magpie.VT.lift                  = VTAero(magpie.VT.disc, lift_VT);
magpie.VT.torque                = VTTorque(magpie.VT.disc, magpie.VT.lift);
[magpie.VT.SF, magpie.VT.BM]    = CalculateSFandBM(magpie.VT.lift, ...
                                    magpie.VT.disc);

save("MagpieStruct","magpie")


plot(magpie.VT.disc, magpie.VT.torque)


%plotting stuff for the repoort

[magpie.HT.lift_landing.SF.n_max, magpie.HT.lift_landing.BM.n_max] = CalculateSFandBM(magpie.HT.lift_landing.n_max, magpie.HT.disc);
[magpie.HT.lift.SF.n_max, magpie.HT.lift.BM.n_max] = CalculateSFandBM(magpie.HT.lift.n_max, magpie.HT.disc);

magpie.HT.inertial.SF.n_max = zeros(1, length(magpie.HT.disc));
magpie.HT.inertial.BM.n_max = zeros(1, length(magpie.HT.disc));

for i = 1:length(magpie.HT.disc)
    [magpie.HT.inertial.BM.n_max(i), magpie.HT.inertial.SF.n_max(i)] = HTweight(HT, magpie.HT.disc(i));
end

magpie.HT.combined_landing.SF.n_max = magpie.HT.lift_landing.SF.n_max - 4.5*magpie.HT.inertial.SF.n_max;
magpie.HT.combined_landing.BM.n_max = magpie.HT.lift_landing.BM.n_max - 4.5*magpie.HT.inertial.BM.n_max;

figure()
plot(magpie.HT.disc, 4.5*magpie.HT.inertial.SF.n_max)
hold on
plot(magpie.HT.disc, magpie.HT.lift.SF.n_max)
hold off
magpie.HT.combined.SF.n_max = magpie.HT.lift.SF.n_max - 3.75*magpie.HT.inertial.SF.n_max;
magpie.HT.combined.BM.n_max = magpie.HT.lift.BM.n_max - 3.75*magpie.HT.inertial.BM.n_max;

figure()
p1 = plot(magpie.HT.disc, magpie.HT.combined.SF.n_max, "-k", "LineWidth", 2, "DisplayName", "Loading at $V_A$ and $V_D$");
hold on
p2 = plot(magpie.HT.disc, magpie.HT.combined_landing.SF.n_max, "-r", "LineWidth", 2, "DisplayName", "Landing nose off");


plot([-1000 1000], [0 0], "-k", "LineWidth", 0.5)
grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")
legend([p1 p2], "Interpreter","latex", "FontSize", 16)
xlabel("Station along HT (m)", "Interpreter","latex", "FontSize", 16)
ylabel("Shear force (N)", "Interpreter","latex", "FontSize", 16)
xlim([0 b_tail/2])

hold off

figure()
p1 = plot(magpie.HT.disc, magpie.HT.combined.BM.n_max, "-k", "LineWidth", 2, "DisplayName", "Loading at $V_A$ and $V_D$");
hold on
p2 = plot(magpie.HT.disc, magpie.HT.combined_landing.BM.n_max, "-r", "LineWidth", 2, "DisplayName", "Landing nose off");


plot([-1000 1000], [0 0], "-k", "LineWidth", 0.5)
grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")
legend([p1 p2], "Interpreter","latex", "FontSize", 16)
xlabel("Station along HT (m)", "Interpreter","latex", "FontSize", 16)
ylabel("Bending Moment (Nm)", "Interpreter","latex", "FontSize", 16)
xlim([0 b_tail/2])

hold off


figure()
p1 = plot(magpie.HT.disc, magpie.HT.torque.n_max, "-k", "LineWidth", 2, "DisplayName", "Loading at $V_A$ and $V_D$");
hold on
p2 = plot(magpie.HT.disc, magpie.HT.torqueLanding.n_max, "-r", "LineWidth", 2, "DisplayName", "Landing nose off");


plot([-1000 1000], [0 0], "-k", "LineWidth", 0.5)
grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")
legend([p1 p2], "Interpreter","latex", "FontSize", 16, "Location", "southeast")
xlabel("Station along HT (m)", "Interpreter","latex", "FontSize", 16)
ylabel("Torque (Nm)", "Interpreter","latex", "FontSize", 16)
xlim([0 b_tail/2])

hold off


figure()
p1 = plot(magpie.VT.disc, magpie.VT.SF, "-k", "LineWidth", 2, "DisplayName", "Loading at $V_A$ and $V_D$");
hold on
plot([-1000 1000], [0 0], "-k", "LineWidth", 0.5)
grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

xlabel("Station along VT (m)", "Interpreter","latex", "FontSize", 16)
ylabel("Shear force (N)", "Interpreter","latex", "FontSize", 16)
xlim([0 b_vtail/2])

hold off

figure()
p1 = plot(magpie.VT.disc, magpie.VT.BM, "-k", "LineWidth", 2, "DisplayName", "Loading at $V_A$ and $V_D$");
hold on
plot([-1000 1000], [0 0], "-k", "LineWidth", 0.5)
grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")
xlabel("Station along VT (m)", "Interpreter","latex", "FontSize", 16)
ylabel("Bending moment (Nm)", "Interpreter","latex", "FontSize", 16)
xlim([0 b_vtail/2])

hold off

figure()
p1 = plot(magpie.VT.disc, magpie.VT.torque, "-k", "LineWidth", 2, "DisplayName", "Loading at $V_A$ and $V_D$");
hold on
plot([-1000 1000], [0 0], "-k", "LineWidth", 0.5)
grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")
xlabel("Station along VT (m)", "Interpreter","latex", "FontSize", 16)
ylabel("Torque (Nm)", "Interpreter","latex", "FontSize", 16)
xlim([0 b_vtail/2])

hold off
