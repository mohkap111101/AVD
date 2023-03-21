%Changes req'd:

% - Think about optimising for rib spacing
% - Investigate t_eff

% clc
% clear all
% close all
load('liftSurfGeom.mat');
load('effectiveThicknesses.mat');
t_e = eff_thicknesses(1:end-1)';
% Now, it is possible to discretise the wing span into increments of 's':
span = 12.165; % Wing Span/2 (m)

%RUN WINGCOVER.MAT
rib_coordinates = rib.distL(1,:)
rib_spacing = rib.distL(2,:)
s=rib_spacing;

s_r=rib_coordinates
% The chord distribution along the span is:

root_chord=4.073949602;
tip_chord=1.756951274;

c_r = ((tip_chord-root_chord)/span).*rib_coordinates + root_chord;
c_r = c_r*0.4;

no_ribs= length(rib_coordinates);

%% RC's Remerz function to find BM at each rib:
% Bending Moments:
BM_i=zeros(1,no_ribs);
BM_a=zeros(1,no_ribs);

for i=1:no_ribs
    location=s_r(i);
    BM_i(i)= Remerz(location);
    BM_a(i)= aeroMom(wing,location);
end

BM=3.75*(BM_a+BM_i);

% Implementing Remi's script to find the coordinates of the wingbox, for a
% chord of 1
figure
coords=wingbox(wing);

figure
hold on
plot(s_r,c_r,'x')
ylabel("Rib chord")
xlabel("Rib span (location)")
legend('Ribs, evenly distributed','Constraints')
hold off

% Now, it is possible to find the wing box height 
h_c=zeros(1,no_ribs); % Depth of wing box along chord (m)

for j=1:length(c_r)
    h_c(j) = c_r(j)*((coords(1,2)-coords(2,2))+(coords(4,2)-coords(3,2)))/2;
end

figure
plot(s_r,h_c,'x')
ylabel("Span")
xlabel("Wing box height")

%% Crushing 

% To determine the rib thickness, it is necessary to equate the yield
% stress of the material to the critical buckling stress of the material

% From Zahra, the critical buckling stress is given by: 3.62E(t/h)^2 where
% E is Young's Modulus, t is rib thickness, and h is wing box depth

% From Daqing, the critical yield stress is given by: F/t*c where F is the
% crush force (see video 11:40), and c is chord at rib position.



% t_e = []; % Effective length (m)
% for l=1:length(c_r)
%     t_e(l)=0.002;
% end
% 
E = 7.35*10^10; % Young's Modulus (Pa) - Used Al currently
I=[];
t_r=[];
% 
% % Crushing force:
% 
I = (c_r .* t_e.^3)./12 + (c_r .* t_e .* (h_c./2).^2);
F = (BM.^2 .* s .* h_c .* t_e .* c_r) ./ (2 .* E .* I.^2) 
t_r = ((F .* h_c.^2)./(3.62.*E.*c_r)).^(1/3);
% 
% figure
plot(s_r,t_r,'x')
% xlabel("Rib Span")
% ylabel("Rib Thickness")

figure()
plot(s_r,I, "-ok", "LineWidth", 1.5)
grid on 
grid minor
xlabel("Rib Position (m)","Interpreter","latex", "FontSize", 16)
ylabel("Second moment of area, $I$, (m$^4$)","Interpreter","latex", "FontSize", 16)
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure()
plot(s_r,F, "-sk", "LineWidth", 1.5)
grid on 
grid minor
xlabel("Rib Position (m)","Interpreter","latex", "FontSize", 16)
ylabel("Crushing force, $F$, (N)","Interpreter","latex", "FontSize", 16)
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure()
plot(s_r,t_r, "-^k", "LineWidth", 1.5)
grid on 
grid minor
xlabel("Rib Position (m)","Interpreter","latex", "FontSize", 16)
ylabel("Rib thickness, $t_r$, (m)","Interpreter","latex", "FontSize", 16)
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

