% Master script for estimating the aerodynamic parameters and stability
% derivatives of the aircraft 

% AVD Group 19 - Project MAGPIE;

% Lara Alves, Shiven Chandarana, Rohan Chotai, Mohammad Kapadia, Remi
% Obasa, Mostafa Trihan 

%% Housekeeping 

clear
clc
close all

%% Constants and parameters 

% Mission Parameters
rho0 = isa_func("rho",0); %density at sea level
R0 = 360*100; % cruise 0 range (first cruise for which values are given) 
h0 = 25000; % height for cruise 0 (ft)
R = [710,420,370].*1000; % Ranges for different cruise legs (m)
C = 1/7200; % Fuel Consumption#
t_turnaround = 30;
TW_idle = 0.3;
LD_max = 19.14; 
M_c0 = 0.69; %Cruise 0 Mach number
V_c0 = M_c0 * isa_func("a",h0); % cruise 0 velocity
V_md_c0 = V_c0/(3^(1/4)); % assuming that v_c1 = v_md
V_imd = V_md_c0*sqrt(isa_func("rho",25000)/rho0); % working out v_imd
MTOW = 45119;

h1 = 25000;
Mc1 = 0.69;
h2 = 28000;
Mc2 = 0.69;
V_c2 = Mc2 * isa_func("a",h2);

% Mass Fractions
w_engine_start = 0.99;
w_taxi = 0.99;
w_takeoff = 0.995;
w_descent = 0.99;
w_landing = 0.992;
w_turnaround = exp(-1*t_turnaround * 60 * C * TW_idle / 9.81);

% Total masses
M0 = 45119;
M_start_c2 = 0.8807*M0;
M_end_c2 = 0.85108156 * M0;
W_cruise_2 = 42980.21236 * 9.81;

% Design point
WS = 5419.78;
TW = 0.277;

% Design Parameters
AR = 9;
e = 0.8;
Cl_max_landing = 2.9;
Cl_max_takeoff = 0.8*Cl_max_landing;
Cl_max = 1.5;
Sref = 85.08;
b = 27.67;
FuselageDiameter = 3;
Ne = 2;
x_cg = 11.65;

% Fuselage Parameters
de_fuselage = 3.34;
D = de_fuselage;
l_fuselage = 27.86;
f_fuselage = l_fuselage/de_fuselage;
L1 = 3.3496;    
L2 = 20.488;
L3 = 2.0183;
Swet_fuselage = pi*D/4 * ( 1/(3*L1^2) * ((4*L1^2 + D^2/4)^1.5 - D^3/8) - D + 4*L2 + 2*sqrt(L3^2 + D^2/4));
A_max_fuselage = (de_fuselage/2)^2 * pi;

eta_fuselage = (de_fuselage/2)/(b/2);
upsweep_angle = 13*pi/180;

% Wing parameters
alpha_stall_wing = 15;
alpha_min_Cl = -10;
x_c_m_wing = 0.30;
t_c_m_wing = 0.12;
quarter_cord_sweep_wing = 0;
max_sweep_wing = 5;
trailing_edge_sweep_wing = 5;
leading_edge_sweep_wing = 5;
c0 = 4.3923;
lambda = 0.4;
c_tip = lambda*c0;
c_bar_wing = 2/3 * (c0^2 * (b/2))/(Sref) * (lambda^2 + lambda + 1); % mean aerodynamic chord
y_mac = 1/3 * (c0 * (b/2)^2)/(Sref) * (2*lambda + 1);
x_cg_wing_qc = 0.4894;
x_wing = 11.5;
x_c_4 = x_wing - x_cg_wing_qc;
if(Mc2 > 0.4)
    x_ac_wing = x_c_4 + 0.26*(Mc2 - 0.4)^2.5 * sqrt(Sref);
else
    x_ac_wing = x_c_4;
end
S_exposed = (b-D)*(c0 + c_tip)/2;
Swet_wing = 1.07*S_exposed*2;
dihedral = 7*pi/180; % Radians - dihedral

% Flaps Parameters
sweep_hinge_line = 5;
Rf = 0.25;
S_flapped = 0.623;
c_dash_c_50_degs = 1.03;
c_dassh_c_30_degs = 1.03;
delta_flap_takeoff = 30;
delta_flap_landing = 50;

% Slats Parameters
sweep_hinge_line_slats = 5;
Rf_slats = 0.17;
c_dash_c_slats = 1.17;
S_slatted = 0.828;
delta_slat = 25.6; 

% Tail Parameters
x_c_m_tail = 0.30;
t_c_m_tail = 0.12;
b_tail = 9.5377;
quarter_cord_sweep_tail = 26;
max_sweep_tail = 35;
trailing_edge_sweep_tail = 35;
leading_edge_sweep_tail = 35;
S_tail = 20.215;
c0_tail = 2.8666666;
lambda_tail = 0.4;
ctip_tail = lambda_tail * c0_tail;
c_bar_tail = 2/3 * (c0_tail^2 * (b_tail/2))/(S_tail) * (lambda_tail^2 + lambda_tail + 1);
y_mac_tail = 1/3 * (c0_tail * (b_tail/2)^2)/(S_tail) * (2*lambda_tail + 1);
AR_tail = 4.5;
h_H = 2.643500128;
eta_H = 0.9;
S_exposed_tail = (c0_tail + ctip_tail)/2 * b_tail;
Swet_tail = 1.05*S_exposed_tail*2;
x_tail = 25.054;
x_cg_tail_le = 1.22;
x_c_4_tail = 24.645; %x_tail - (x_cg_tail_le - 0.25*c0_tail);
if(Mc2 > 0.4)
    x_ac_h = x_c_4_tail  + 0.26*(Mc2 - 0.4)^2.5 * sqrt(S_tail);
else
    x_ac_h = x_c_4_tail;
end
l_H = x_ac_h - x_ac_wing;

% Vertical Tail Parameters

x_c_m_vtail = 0.30;
t_c_m_vtail = 0.12;
b_vtail = 5.5848*2;
quarter_cord_sweep_vtail = 35;
mid_chord_sweep_vtail = 33.7;
max_sweep_vtail = 45;
trailing_edge_sweep_vtail = 40;
leading_edge_sweep_vtail = 40;
S_vtail = 17.328;
c0_vtail = 4.153;
lambda_vtail = 0.4;
ctip_vtail = lambda_vtail * c0_vtail;
c_bar_vtail = 2/3 * (c0_vtail^2 * (b_vtail/2))/(S_vtail*2) * (lambda_vtail^2 + lambda_vtail + 1);
y_mac_vtail = 1/3 * (c0_vtail * (b_vtail/2)^2)/(S_vtail*2) * (2*lambda_vtail + 1);
AR_vtail = 3.6;
S_exposed_vtail = (c0_vtail + ctip_vtail)/2 * b_vtail/2;
Swet_vtail = 1.05*S_exposed_vtail*2;
x_c_4_vtail = 23.253;
if(Mc2 > 0.4)
    x_ac_vtail = 23.3543 + 0.26*(Mc2 - 0.4)^2.5 * sqrt(S_vtail*2);
else
    x_ac_vtail = x_c_4_vtail;
end
x_vtail = 26.7957;
l_vt_cg = x_ac_vtail - x_cg;
Z_vt = 4.17175;

% Nacelle Parameters
de_nacelle = 1.5;
l_nacelle = 3;
f_nacelle = l_nacelle/de_nacelle;
x_nacelle_fuselage = 1.5;
Swet_nacelle = pi * (de_nacelle/2) * l_nacelle;

% Landing gear parameters - Nose
d = 0.6477;
w = 0.2;
Ref_Area = d*w;
e = 1.06;
a_nose_to_gear = 2.7;
a_d = a_nose_to_gear/d;
del_Cds_Nose_gear = 0.5;

% Main landing gear parameters
d_wheel = 0.755;
width_wheel = 0.29;
d_mg = 0.635 + 0.755/2;
w_mg = 0.2184 + 0.29;
ref_area_mg = d_mg * w_mg;
S_frontal = (0.635*0.2184) + (0.755*0.29);


% Values at different altitudes
[T_leg1, a_leg1, P_leg1, rho_leg1] = atmosisa(distdim(h1, "ft", "m"));
[T_leg2, a_leg2, P_leg2, rho_leg2] = atmosisa(distdim(h2, "ft", "m"));


% Calculating the Reynolds number
Re_crit = 10^5;
Re_wing = ReynoldsNumber(rho_leg2, a_leg2, T_leg2, c_bar_wing, Mc2);
Re_fuselage = ReynoldsNumber(rho_leg2, a_leg2, T_leg2, de_fuselage, Mc2);
Re_nacelle = ReynoldsNumber(rho_leg2, a_leg2, T_leg2, de_nacelle, Mc2);
Re_tail = ReynoldsNumber(rho_leg2, a_leg2, T_leg2, c_bar_tail, Mc2);
Re_vtail = ReynoldsNumber(rho_leg2, a_leg2, T_leg2, c_bar_vtail, Mc2);

% Reynolds numbers for different chords for wing
Re_wing_root = ReynoldsNumber(rho_leg2, a_leg2, T_leg2, c0, Mc2);
Re_wing_tip = ReynoldsNumber(rho_leg2, a_leg2, T_leg2, c_tip, Mc2);

% Reynolds numbers for different chords for tail

Re_tail_root = ReynoldsNumber(rho_leg2, a_leg2, T_leg2, c0_tail, Mc2);
Re_tail_tip = ReynoldsNumber(rho_leg2, a_leg2, T_leg2, ctip_tail, Mc2);

% Reynolds numbers for different chords for vertical stabiliser (vtail)#

Re_vtail_root = ReynoldsNumber(rho_leg2, a_leg2, T_leg2, c0_vtail, Mc2);
Re_vtail_tip = ReynoldsNumber(rho_leg2, a_leg2, T_leg2, ctip_vtail, Mc2);

Cl_cruise = M_start_c2*9.81/(0.5*rho_leg2*(Mc2*a_leg2)^2*Sref);

M_dd = 0.87/cosd(quarter_cord_sweep_wing) - t_c_m_wing/(cosd(quarter_cord_sweep_wing))^2 - ((Cl_cruise/0.9)/10)/(cosd(quarter_cord_sweep_wing))^3;

%% Thrust lapse model
% TR =  1.072;
% 
% A = 0:10000:40000;
% Mach = 0:0.001:1;
% 
% [Altitudes,M] = meshgrid(A, Mach);
% 
% for j=1:length(A)
%     for i=1:length(Mach)
%         ThrustRatios(j,i) = ThrustLapse(A(j), TR, Mach(i), 1.4);
%     end
% end
% 
% figure()
% for j = 1:length(A)
%     plot(Mach, ThrustRatios(j,:), "LineWidth", 2)
%     hold on
% end
% 
% legend("0 ft - Sea level", "10000 ft", "20000 ft", "30000 ft", "40000 ft", "Interpreter","latex", "FontSize", 16)
% xlabel("Mach number", "Interpreter","latex", "FontSize", 16)
% ylabel("Thrust Ratio", "Interpreter","latex", "FontSize", 16)
% set(gca,'FontSize', 16, "TickLabelInterpreter", "Latex")

%% Plot surface for thrust lapse model 
% surface = surf(A, Mach, ThrustRatios');
% surface.EdgeColor = 'none';
% 
% xlabel("Altitude", "Interpreter","latex", "FontSize", 16)
% ylabel("Mach number", "Interpreter","latex", "FontSize", 16)
% zlabel("Thrust Ratio", "Interpreter","latex", "FontSize", 16)
% 
% set(gca,'FontSize', 16, "TickLabelInterpreter", "Latex")

%% Lift Estimation

% Wing lift estimation 

M_crit = M_dd - (0.1/80)^(1/3);

% Calculate the lift curve slope of the wing
a_theory = 2*pi + 4.7*t_c_m_wing*(1 + 0.00375*-5);

beta = sqrt(1-Mc2^2);
a_atheory = 0.95;
a_sectional = (1.05/beta) * (a_atheory) * a_theory;

eta_airfoil = 0.92;%a_sectional * beta/(2*pi);

F_d_b_factor = (S_exposed/Sref) * 1.07*(1 + D/b)^2;

if(F_d_b_factor >= 1)
    F_d_b_factor = 0.98;
end

if(Mc2 > M_crit)
    a_wing = ((2*pi*AR)/(2 + sqrt((AR^2 * beta^2)/(eta_airfoil^2) * (1 + ((tand(quarter_cord_sweep_wing))^2)/(beta^2)) + 4))) * F_d_b_factor;
    a_wing_M0  = ((2*pi*AR)/(2 + sqrt((AR^2 * 1^2)/(eta_airfoil^2) * (1 + ((tand(quarter_cord_sweep_wing))^2)/(1^2)) + 4))) * F_d_b_factor;
else
    beta = beta;
    a_wing = ((2*pi*AR)/(2 + sqrt((AR^2 * beta^2)/(eta_airfoil^2) * (1 + ((tand(quarter_cord_sweep_wing))^2)/(beta^2)) + 4))) * F_d_b_factor;
    a_wing_M0  = ((2*pi*AR)/(2 + sqrt((AR^2 * 1^2)/(eta_airfoil^2) * (1 + ((tand(quarter_cord_sweep_wing))^2)/(1^2)) + 4))) * F_d_b_factor;
end


  % Zero lift angle of attack - same as the aeorfoil used for the wing and
  % tail
alpha_0_wing = -3; % degrees
alpha_0_tail = 0;

  % Calculating CLmax 

% Using potential flow solver to calculate the different Cl distribution
% at different angles of attack. Then finding the angle at which the Clmax
% distribution touches the Cl distribution and then integrating that curve over
% the span of the wing.

Clmax_aerofoil = 1.65;
AoAs = [-2, 0, 2, 4, 6, 8, 10, 12, 14, 16];
Cls = [0.1, 0.3, 0.55, 0.75, 0.95, 1.15, 1.35, 1.5, 1.6, 1.65];

figure()

for i = 1:length(Cls)
    etas = linspace(0, 1, 500);
    cs = chord_dist(etas, lambda, c0);
    CLs = SectionaCL(Cls(i), cs, b, Sref);
    
    plot(etas, CLs)
    hold on
end

plot([0, 1], [Clmax_aerofoil, Clmax_aerofoil])
plot([eta_fuselage eta_fuselage], [0, 20])

ylim([0 2.5])

% Tail Lift estimation

a_theory_tail = 2*pi + 4.7*t_c_m_wing*(1 + 0.00375*-5);

beta = sqrt(1-Mc2^2);
a_atheory_tail = 0.95;
a_sectional_tail = (1.05/beta) * (a_atheory) * a_theory;

eta_airfoil_tail = 0.92;%a_sectional_tail * beta/(2*pi);

if(Mc2>M_crit)
    a_tail = ((2*pi*AR_tail)/(2 + sqrt((AR_tail^2 * beta^2)/(eta_airfoil_tail^2) * (1 + ((tand(max_sweep_tail))^2)/(beta^2)) + 4)))*F_d_b_factor;
else
    beta = beta;
    a_tail = ((2*pi*AR_tail)/(2 + sqrt((AR_tail^2 * beta^2)/(eta_airfoil_tail^2) * (1 + ((tand(max_sweep_tail))^2)/(beta^2)) + 4)))*F_d_b_factor;
end

% Lift curve slope of the vertical stabiliser

a_theory_tail = 2*pi + 4.7*t_c_m_wing*(1 + 0.00375*-5);

beta = sqrt(1-Mc2^2);
a_atheory_vtail = 0.95;
a_sectional_vtail = (1.05/beta) * (a_atheory) * a_theory;

eta_airfoil_vtail = 1;%a_sectional_tail * beta/(2*pi);

if(Mc2>M_crit)
    a_vtail = ((2*pi*AR_vtail)/(2 + sqrt((AR_vtail^2 * beta^2)/(eta_airfoil_vtail^2) * (1 + ((tand(mid_chord_sweep_vtail))^2)/(beta^2)) + 4)));
else
    beta = beta;
    a_vtail = ((2*pi*AR_vtail)/(2 + sqrt((AR_vtail^2 * beta^2)/(eta_airfoil_vtail^2) * (1 + ((tand(mid_chord_sweep_vtail))^2)/(beta^2)) + 4)));
end

% Whole aircraft lift curve slope estimation

KA = 1/AR - 1/(1 + AR^1.7);
K_lambda = (10 - 3*lambda)/7;
Kh = (1-abs(h_H/b))/((2*l_H/b)^(1/3));

depsilon_dalpha = 4.44*(KA * K_lambda * Kh * sqrt(cos(0)))^1.19 * a_wing/a_wing_M0;

Cl_max_aircraft = a_wing*alpha_stall_wing*(pi/180) + S_tail/Sref * eta_H * (1-depsilon_dalpha)*a_tail*alpha_stall_wing*(pi/180);
Cl_min_aircraft = a_wing*alpha_min_Cl*(pi/180) + S_tail/Sref * eta_H * (1-depsilon_dalpha)*a_tail*alpha_min_Cl*(pi/180);

%% Lift estimation - Takeoff and landing 
% Raymer

del_alpha_zero_takeoff_aerofoil = -10; % Change in zero lift angle of attack in takeoff - degrees
del_alpha_zero_landing_aerofoil = -15; % Change in zero lift angle of attach in landing - degrees

del_alpha_zero_takeoff = del_alpha_zero_takeoff_aerofoil * S_flapped* cosd(sweep_hinge_line);
del_alpha_zero_landing = del_alpha_zero_landing_aerofoil * S_flapped* cosd(sweep_hinge_line);

% Estimating lift curve slope increase - landing
%a_slats_extended = a_wing * (1 + (c_dash_c_slats - 1)*(S_slatted));
a_flaps_extended_takeoff = a_wing * (1 + (c_dassh_c_30_degs - 1)*(S_flapped));
a_flaps_extended_landing = a_wing * (1 + (c_dash_c_50_degs - 1)*(S_flapped));

del_CLmax_landing_flaps = 0.92;
del_CLmax_landing_slats = 0.643964;

del_CLmax_takeoff_flaps = del_CLmax_landing_flaps*0.85*0.6;

a_wing = a_wing;
alpha_0_wing = alpha_0_wing;% + del_alpha_zero_takeoff;

a_total = a_wing + (S_tail/Sref)*a_tail*(1-depsilon_dalpha);

a_tail = (S_tail/Sref)*a_tail*(1-depsilon_dalpha);
Lift_tail = a_tail * (2.2 - 0.25) * pi/180 * 0.5 * rho_leg2 * V_c2^2 * S_tail;

alpha_c = ((1/(a_total * pi/180)) * (9.81*(M_start_c2))/(0.5*rho_leg2 * (Mc2*a_leg2)^2 * Sref)) + alpha_0_wing;

% Estimating the relative incidences of the wing and fuselage

phi_g = -2; % Geometic washout of the wing 
delta_phi_mgc = ((1 + 2*lambda)/(3 + 3*lambda)) * phi_g;

alpha_FOPT = 2;

i_w = 1.0;

%% Drag Estimation - Cruise

% List indexes - [1 - wing, 2 - fuselage, 3 - nacelle, 4 - H tail, 5 - V Tail]
Swet = [Swet_wing, Swet_fuselage, Swet_nacelle, Swet_tail, Swet_vtail];
FF = [];
Q = [];

% Calculate Cf from Re for each component
% Cf(1) = skinFrictionCoefficient(Re_wing, Re_crit, Mc2);
% Cf(2) = skinFrictionCoefficient(Re_fuselage, Re_crit, Mc2);
% Cf(3) = 2*skinFrictionCoefficient(Re_nacelle, Re_crit, Mc2);
% Cf(4) = skinFrictionCoefficient(Re_tail, Re_crit, Mc2);
% Cf(5) = skinFrictionCoefficient(Re_vtail, Re_crit, Mc2);

% Estimate the form drag for each component
FF(1) = FormFactor_WingSections(t_c_m_wing, x_c_m_wing, quarter_cord_sweep_wing, Mc2);
FF(2) = FormFactorOtherComps(f_fuselage, "fuselage");
FF(3) = FormFactorOtherComps(f_nacelle, "nacelle");
FF(4) = FormFactor_WingSections(t_c_m_tail, x_c_m_tail, quarter_cord_sweep_tail, Mc2);
FF(5) = FormFactor_WingSections(t_c_m_vtail, x_c_m_vtail, quarter_cord_sweep_vtail, Mc2);

% Estimate the Interference factor for each compoenent
Q(1) = 1;
Q(2) = 1;
Q(3) = InterferenceDrag(de_nacelle, x_nacelle_fuselage);
Q(4) = 1.05;
Q(5) = 1.03;

% Induced Drag#
%Cdi = (1/pi*AR*e)*Cl^2;


% Drag Estimation by components of drag - Gudmundson

% Skin Friction drag coefficient using the mixed turbulent and laminar BL
% model

X_tr_wing = 0.55;

X_tr_tail = 0.02;

X_tr_vtail = 0.5;

X_tr_fuselage = 0.90; 

X_tr_nacelle = 0.5;

Cf_wing = FindCf(X_tr_wing, Re_wing, Mc2, M_crit);
Cf_tail = FindCf(X_tr_tail, Re_tail, Mc2, M_crit);
Cf_vtail = FindCf(X_tr_vtail, Re_vtail, Mc2, M_crit);
Cf_fuselage = FindCf(X_tr_fuselage, Re_fuselage, Mc2, M_crit);
Cf_nacelle = 2* FindCf(X_tr_nacelle, Re_nacelle, Mc2, M_crit);

Cf = [Cf_wing, Cf_fuselage, Cf_nacelle, Cf_tail, Cf_vtail];

% Calculaiting Cd0

Cd_0_skin_friction = sum(Cf.*FF.*Q.*Swet)/Sref;


% Wave drag estimation 
Cd_maxD = 0.03;
M_maxD = 1.2;

A = (atanh((2*Cd_maxD - 0.0002)/(Cd_maxD) - 1) - atanh((0.0002)/(Cd_maxD) - 1))/(M_maxD - M_crit);
B = atanh((0.0002)/(Cd_maxD) - 1) - A*M_crit;

Cd_wave = Cd_maxD/2 * (1 + tanh(A*Mc2 + B));
if(Mc2 > M_crit)
    Cd_wave_locks = 20*(Mc2-M_crit)^4;
else
    Cd_wave_locks = 0;
end
Cd_upsweep = 3.83*(pi*D^2/(4*Sref))*upsweep_angle^2.5;

CD0_cruise = Cd_0_skin_friction + Cd_wave_locks + Cd_upsweep;

a_delta_e = 6;
delta_e = 0;

% Calculate the oswald efficiency - e
f_wing = f(lambda);
e_wing = 1/((1+0.12*Mc2^6)*(1+0.142 + f_wing*AR*(10*t_c_m_wing)^0.33 + (0.1*(3*Ne + 1))/((4+AR)^0.8)));
s_wing = 1 - 2*(D/b)^2;
e_theo_wing = 1/(1+f_wing*AR);
K = 0.38;

e_wing = 1/((1/e_theo_wing*s_wing)  + K*CD0_cruise*pi*AR);

f_tail = f(lambda_tail);
e_tail = 1/((1+0.12*Mc2^6)*(1+0.142 + f_tail*AR_tail*(10*t_c_m_tail)^0.33 + (0.1*(3*Ne + 1))/((4+AR_tail)^0.8)));
s_tail = 1 - 2*(0.4/b_tail)^2;
e_theo_tail = 1/(1+f_tail*AR_tail);

e_tail = 1/((1/e_theo_tail*s_tail)  + K*CD0_cruise*pi*AR_tail);



%% Drag Estimation - Takeoff and Landing
% 
del_Cd_landing_gear_nose = (d*w)/Sref * del_Cds_Nose_gear;
del_Cd_retracted_gear = 0.04955*exp(5.615*S_frontal/(d_mg*w_mg)) * (d_mg*w_mg/Sref);

% Young's method
del1_function = 179*Rf^4 - 111.6*Rf^3 + 28.929*Rf^2 + 2.3705*Rf - 0.0089;
del2_function_takeoff = -3.9877e-12 * delta_flap_takeoff^6 + 1.1685e-9*delta_flap_takeoff^5 - 1.2846e-7*delta_flap_takeoff^4 + 6.1742e-6*delta_flap_takeoff^3 + -9.89444e-5*delta_flap_takeoff^2 + 6.8324e-4*delta_flap_takeoff - 3.892e-4;
del2_function_landing = -3.9877e-12*delta_flap_landing^6 + 1.1685e-9*delta_flap_landing^5 - 1.2846e-7*delta_flap_landing^-4 + 6.1742e-6*delta_flap_landing^3 + -9.89444e-5*delta_flap_landing^2 + 6.8324e-4*delta_flap_landing - 3.892e-4;

del_Cd_flaps_takeoff = del1_function * del2_function_takeoff * S_flapped/Sref;
del_Cd_flaps_landing = del1_function*del2_function_landing * S_flapped/Sref;

% Raymer 
F_flap = 0.0074; % From raymer for slotted flaps
del_Cd0_flap_takeoff = F_flap*Rf*S_flapped * (delta_flap_takeoff-10);
del_Cd0_flap_landing = F_flap*Rf*S_flapped * (delta_flap_landing-10);

CD0_cruise = CD0_cruise + del_Cd_retracted_gear + del_Cd_landing_gear_nose + del_Cd0_flap_landing;

% Change in induced drag due to flaps 
kf = 0.28;
del_Cdi_flaps = kf^2*((del_CLmax_landing_flaps + del_CLmax_landing_slats)/1.3^2)^2*cos(0);

%% Stability - Longitudinal

% Figure out the static margin SM of the aircraft using Cm_alpha and
% CL_alpha 

Zt = 2.34;

K_fus = 0.4;
K_nacelle = 0.4;
CM_alpha_fuselage = (K_fus * de_fuselage^2 * l_fuselage)/(c_bar_wing*Sref);
CM_alpha_nacelle = (K_nacelle * de_nacelle^2 * l_nacelle)/(c_bar_wing*Sref);

x_np_c_bar = (a_wing * x_ac_wing/c_bar_wing - CM_alpha_fuselage - 2*CM_alpha_nacelle + eta_H * a_tail*(1-depsilon_dalpha)*S_tail/Sref * x_ac_h/c_bar_wing)/(a_wing + eta_H * a_tail*(1-depsilon_dalpha)*S_tail/Sref);
x_np = x_np_c_bar * c_bar_wing;

SM = x_np_c_bar - x_cg/c_bar_wing;

d_CM_d_alpha = -a_total*SM;

dCM_cg_d_alpha = -a_wing * (x_ac_wing - x_cg)/(c_bar_wing) + CM_alpha_fuselage + 2*CM_alpha_nacelle - eta_H * a_tail*(1-depsilon_dalpha)*(S_tail/Sref)*(x_ac_h - x_cg)/c_bar_wing;

C_m_0_airfoil = -0.07;
C_M_0_wing = (C_m_0_airfoil * ((AR*cos(quarter_cord_sweep_wing)^2)/(AR + 2*cos(quarter_cord_sweep_wing))) - 0.01*phi_g) * (a_wing)/a_wing_M0;

C_M_0_wing_landing = C_M_0_wing - 0.2655;
C_M_0_wing_takeoff = C_M_0_wing - 0.18216;

alphas = 7.6:.01:7.7;
i_hs = 3.5:.01:3.6;

for i = 1:length(i_hs)
    for j = 1:length(alphas)
        
        CLw = a_wing*(deg2rad(alphas(j)) + deg2rad(i_w) - deg2rad(alpha_0_wing));
        CLH = a_tail*((deg2rad(alphas(j)) + deg2rad(i_w) - deg2rad(alpha_0_wing)) * (1-depsilon_dalpha) + (deg2rad(i_hs(i)) - deg2rad(i_w)) - (deg2rad(alpha_0_tail) - deg2rad(alpha_0_wing))) + a_delta_e * deg2rad(delta_e); 
        
        [Cd, Cdi] = GetTotalDragCoeff(CD0_cruise, alphas(j), i_hs(i), a_wing, a_tail, alpha_0_wing, alpha_0_tail, i_w, depsilon_dalpha, a_delta_e, delta_e, AR, AR_tail, e_wing, e_tail, eta_H, S_tail, Sref);
        q = 0.5 * rho_leg2 * Mc2^2*a_leg2^2;
        T = q*Sref*(Cd);

        Cl(i, j) = CLw + eta_H * S_tail/Sref * CLH;

        CMcg(i, j) = -CLw * (x_ac_wing - x_cg)/c_bar_wing + C_M_0_wing_landing + CM_alpha_fuselage * deg2rad(alphas(j)) - eta_H*CLH*S_tail/Sref * (x_ac_h - x_cg)/c_bar_wing + Zt*T/(q*Sref*c_bar_wing);
        
    end
end
figure()

for i = 1:length(i_hs)

    Cmcg = CMcg(i, :);
    CL = Cl(i, :);
    
    plot(CL, Cmcg, "-ok", "LineWidth", 1.5);

    hold on
end


for i = 1:length(alphas)

    Cmcg = CMcg(:, i);
    CL = Cl(:, i);
    
    plot(CL, Cmcg, "--ok", "LineWidth", 1.5);

    hold on
end


CL_design = W_cruise_2/(0.5 * rho_leg2 * Mc2^2*a_leg2^2 * Sref);

plot([CL_design, CL_design], [-10, 10], "-k", "LineWidth", 1.5);
plot([-10, 10], [0, 0], "-k", "LineWidth", 1.5)

xlim([0.2, 0.6])
ylim([-0.04, 0.05])

xlabel("Lift Coefficient", "Interpreter","latex", "FontSize", 16)
ylabel("Pitching Moment Coefficient", "Interpreter","latex", "FontSize", 16)

set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")


% Cd = GetTotalDrag(CD0_cruise, alpha, i_h)

%% Drag polar of the aircraft and CD and CL vs aoa at cruise 

i_h = 3.98;

aoa_s = -10:0.1:20;
CLw = a_wing*(deg2rad(aoa_s) + deg2rad(i_w) - deg2rad(alpha_0_wing));
CLH = a_tail*((deg2rad(aoa_s) + deg2rad(i_w) - deg2rad(alpha_0_wing)) * (1-depsilon_dalpha) + (deg2rad(i_h) - deg2rad(i_w)) - (deg2rad(alpha_0_tail) - deg2rad(alpha_0_wing))) + a_delta_e * deg2rad(delta_e);

CL_s = CLw + eta_H * S_tail/Sref * CLH;

CD0 = GetParasiticDrag(Mc2, h2, M_dd);
CD_s = GetTotalDragCoeff(CD0, aoa_s, i_h, a_wing, a_tail, alpha_0_wing, alpha_0_tail, i_w, depsilon_dalpha, a_delta_e, delta_e, AR, AR_tail, e_wing, e_tail, eta_H, S_tail, Sref);

figure()
plot(aoa_s, CL_s, "-k", "LineWidth",1.5)
hold on 
plot([-10 20], [0.9*1.65, 0.9*1.65])
xlabel("Angle of Attack, $\alpha$", "Interpreter","latex", "FontSize", 16)
ylabel("Lift Coefficient, $C_L$", "Interpreter","latex", "FontSize", 16)

grid on
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure()
plot(aoa_s, CD_s, "-k", "LineWidth",1.5)
hold on 
plot([-10 20], [0 0], "-k")
plot([0 0], [-0.1 0.6], "-k")
xlabel("Angle of Attack, $\alpha$", "Interpreter","latex", "FontSize", 16)
ylabel("Drag Coefficient, $C_D$", "Interpreter","latex", "FontSize", 16)

grid on
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure()
plot(CD_s, CL_s, "-k", "LineWidth",1.5)
hold on
plot([-0.1 0.3], [0 0], "-k")
plot([0 0], [-1 3], "-k")
xlabel("Drag Coefficient, $C_d$", "Interpreter","latex", "FontSize", 16)
ylabel("Lift Coefficient, $C_L$", "Interpreter","latex", "FontSize", 16)

grid on
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure()
plot(aoa_s, CL_s./CD_s, "-k", "LineWidth",1.5)
hold on
plot([-10 20], [0 0], "-k")
plot([0 0], [-25 25], "-k")
xlabel("Angle of Attack, $\alpha$", "Interpreter","latex", "FontSize", 16)
ylabel("Lift to Drag, $L/D$", "Interpreter","latex", "FontSize", 16)

grid on
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

L_D_max_estimated = max(CL_s./CD_s);
%% Static Stability - Lateral - Raymer taken from USAF Datcom

% Z_wf = -0.2;
% Df = l_fuselage;
% Wf = D;
% 
% 
% CL = CL_design;
% C_n_beta_wing = CL^2 * (1   /(4*pi*AR) - (tan(quarter_cord_sweep_wing)/(pi*AR*(AR + 4*cos(quarter_cord_sweep_wing))) ) * (cos(quarter_cord_sweep_wing) - AR/2 - AR^2/(8*cos(quarter_cord_sweep_wing)) + (6*(x_ac_wing/c_bar_wing - x_cg/c_bar_wing)*sin(quarter_cord_sweep_wing))/(AR) ));
% 
% C_l_beta = 0.0115*dihedral;
% 
% Cl_beta_wing_CL = -0.02; % estimation from graph so could be a lot better 
% 
% C_l_beta_dihedral = -a_wing*dihedral/4 * ((2*(1+2*lambda))/(3*(1+lambda)));
% 
% Cl_beta_wf = -1.2 * (sqrt(AR) * Z_wf * (Df + Wf))/(b^2);
% 
% C_l_beta_wing = (Cl_beta_wing_CL) * CL + C_l_beta_dihedral + Cl_beta_wf;


%% Lateral and Directional Stability using Gudmundsson

% Directional Stability - C_N_beta > 0
CL_s = 3.03/1.1^2;
eta_vtail = 0.95;
B = 15*pi/180; % side slip angle 
C_N_beta_wing = CL_s^2 * (1/(4*pi*AR));

C_N_beta_VT = S_vtail * l_vt_cg/(Sref*b) * a_vtail * B * eta_vtail;

k2_k1_fuselage = 1 - 1/f_fuselage + ((0.24*f_fuselage^3 - 5.6*f_fuselage^2 + 44*f_fuselage - 72)/(1000));

int_fus = 171.9499619;

C_N_beta_fus = -(pi* k2_k1_fuselage * int_fus)/(2*Sref*b);

k2_k1_nacelle = 1 - 1/f_nacelle + ((0.24*f_nacelle^3 - 5.6*f_nacelle^2 + 44*f_nacelle - 72)/(1000));

int_nac = l_nacelle * de_nacelle^2 * pi^2/16;

C_N_beta_nac = -(pi* k2_k1_nacelle * int_nac)/(2*Sref*b);

C_N_beta = C_N_beta_wing + C_N_beta_VT + C_N_beta_nac + C_N_beta_fus;

% Lateral stability - C_L_beta < 0

kappa = 1;

C_L_beta_wing_gamma_zero = CL_s * (-kappa/AR * ((0.71*lambda + 0.29)/(1 + lambda)) + 0.05);
C_L_beta_wing_zw = 0.045841;
C_L_beta_wing_gamma = (((-2.786 - 50.46*AR + 2.653*AR^2) * 10^(-6) ) * dihedral*180/pi ) * pi/180;

C_L_beta_tail_gamma_zero = CL_s * (-kappa/AR_tail * ((0.71*lambda_tail + 0.29)/(1 + lambda_tail)) + 0.05);
C_L_beta_tail_zw = 0;
C_L_beta_tail_gamma = (((-2.706 - 47.47*AR_tail + 4.618*AR_tail^2) * 10^(-6) ) * 0 ) *pi/180;

C_L_beta_VT = -C_N_beta_VT * Z_vt/l_vt_cg + -0.001967;

C_L_beta_wing = C_L_beta_wing_gamma_zero + C_L_beta_wing_zw + C_L_beta_wing_gamma;
C_L_beta_tail = C_L_beta_tail_gamma_zero + C_L_beta_tail_zw + C_L_beta_tail_gamma;

C_L_beta = (C_L_beta_wing) + (C_L_beta_tail) + C_L_beta_VT;

%% Save workspace

save("Mohkap_workspace.mat")

%% Plotting static margins against weights

SM = [0.118918753
0.094328124
0.093369752
0.092482448
0.091989459
0.114628697
0.114944329
0.090333728
0.089293838
0.087479655
0.086936133
0.109534618
0.109747399
0.085110683
0.0839642
0.082679142
0.082087129
0.10464651
0.104760599
0.05184729
0.050598524
0.04930131
0.048657517
0.047126267
0.046460504
0.097211942];


Weights = [44221.2802
44000.1738
43120.17032
42336.23423
41912.87189
41577.56891
40754.7191
40550.94551
39739.9266
38400.0777
38016.07693
37711.94831
36965.6019
36780.77389
36045.15841
35254.83652
34902.28815
34623.06985
33937.85455
33768.16528
33092.80197
32419.26162
32095.06901
31349.42842
31035.93414
30787.64667];

plot(Weights(1:6), SM(1:6), "-o", "Color", "red", "LineWidth", 1.5, "DisplayName","Cruise 1");
hold on
plot(Weights(6:12), SM(6:12), "-o", "Color", "blue", "LineWidth", 1.5, "DisplayName","Cruise 2");
plot(Weights(12:18), SM(12:18), "-o", "Color", "green", "LineWidth", 1.5, "DisplayName","Cruise 3");
plot(Weights(18:23), SM(18:23), "-o", "Color", "black", "LineWidth", 1.5, "DisplayName","Diversion Cruise");
plot(Weights(23:end), SM(23:end), "--o", "Color", "black", "LineWidth", 1.5, "DisplayName","Loiter and Land");
set(gca, "XDir", "reverse")

xlabel("Weight, $W$ (kg)", "Interpreter","latex", "FontSize", 16)
ylabel("Power Static Margin, $SM_{power\,off}$", "Interpreter","latex", "FontSize", 16)

grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")
legend("Interpreter","latex", "FontSize", 16)
