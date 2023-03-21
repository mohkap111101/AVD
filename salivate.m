% Housekeeping
clear
clc

% Load structs
load('tail.mat')
load('liftSurfGeom.mat','HT')
geometry = HT;
upper = HTskin;
rib = HTrib;
rib.distL = [rib.distL, [geometry.span/2-geometry.TEfus; 0]];
distL   = rib.distL(:,[1,3,5,7,end]);
upper.distT(2,upper.distT(1,:)>1) = 0.001;
% p = polyfit(upper.distT(1,:),upper.distT(2,:),6);
% upper.distT(2,:) = polyval(p,upper.distT(1,:));
uT      = upper.distT(2,:);
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

save('MoKapGetsPegged.mat','upper','uT')