% Housekeeping
clear
clc

% Load structs
load('addStraightLines.mat')
load('liftSurfGeom.mat','wing')
geometry = wing;
rib.distL = [rib.distL, [geometry.span/2-geometry.TEfus; 0]];
distL   = rib.distL(:,[1,5,8,12,15,19,22,24,end]);
p = polyfit(upper.distT(1,:),upper.distT(2,:),6);
upper.distT(2,:) = polyval(p,upper.distT(1,:));
p = polyfit(lower.distT(1,:),lower.distT(2,:),6);
lower.distT(2,:) = polyval(p,lower.distT(1,:));
uT      = upper.distT(2,:);
lT      = lower.distT(2,:);
for i = 1:length(distL)-1
    i1          = (upper.distT(1,:) >= distL(1,i));
    i2          = (upper.distT(1,:) <= distL(1,i+1));
    index       = i1 == i2;
    uT(index)   = max(upper.distT(2,index));
    lT(index)   = max(lower.distT(2,index));
end

% Skin process
% Upper
figure()
plot(upper.distT(1,:),upper.distT(2,:), "-.k", "LineWidth", 1.5, DisplayName="Required thicknesss")
hold on
plot(upper.distT(1,:),uT, "-k", "LineWidth", 2, DisplayName="Set Thickness")
xlabel("Span station (m)", Interpreter="latex", FontSize=16)
ylabel("Upper Skin Thickness (m)", Interpreter="latex", FontSize=16)

legend("Interpreter","latex", "FontSize", 16)
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

grid on 
grid minor

% Lower
figure()
plot(lower.distT(1,:),lower.distT(2,:), "-.k", "LineWidth", 1.5, DisplayName="Required thicknesss")
hold on
plot(lower.distT(1,:),lT, "-k", "LineWidth", 2, DisplayName="Set Thickness")
xlabel("Span station (m)", Interpreter="latex", FontSize=16)
ylabel("Lower Skin Thickness (m)", Interpreter="latex", FontSize=16)
legend("Interpreter","latex", "FontSize", 16)
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

grid on 
grid minor

% Rib spacing
figure()
plot(rib.distL(1,:),rib.distL(2,:))

save('uuuugh.mat','upper','lower','uT','lT')