% Housekeeping
clear
clc

% Load structs
load('tail.mat')
rib = HTrib;
upper = HTskin;

distL   = rib.distL(:,[1,3,5, end]);

uT      = upper.distT(2,:);
upper.distT(2,upper.distT(1,:)>1) = 0.001;
% lT      = lower.distT(2,:);
for i = 1:length(distL)-1
    i1          = (upper.distT(1,:) >= distL(1,i));
    i2          = (upper.distT(1,:) <= distL(1,i+1));
    index       = i1 == i2;
    uT(index)   = max(upper.distT(2,index));
%     lT(index)   = max(lower.distT(2,index));
end

% Skin process
% Upper
figure()
plot(upper.distT(1,:),upper.distT(2,:))
hold on
plot(upper.distT(1,:),uT)

% % Lower
% figure()
% plot(lower.distT(1,:),lower.distT(2,:))
% hold on
% plot(lower.distT(1,:),lT)

% Rib spacing
figure()
plot(rib.distL(1,:),rib.distL(2,:))

figure()
p1 = plot(upper.distT(1,:),upper.distT(2,:), "-.k", "LineWidth", 1.5, DisplayName="Required thicknesss");
hold on
p2 = plot(upper.distT(1,:),uT, "-k", "LineWidth", 2, DisplayName="Set Thickness");
xlabel("Span station (m)", Interpreter="latex", FontSize=16)
ylabel("Skin Thickness (m)", Interpreter="latex", FontSize=16)
plot([-1000 1000], [1e-3 1e-3], "-k", LineWidth=0.5)
legend([p1 p2], "Interpreter","latex", "FontSize", 16)
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

grid on 
grid minor

xlim([0 5])

ylim([0.8e-3 1.4e-3])