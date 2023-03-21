clear
clc
close all

load wingCover.mat
load('addStraightLines.mat')


distL   = rib.distL(:,[1,5,8,12,15,19,22,24,end]);
uT      = upper.distT(2,:);
lT      = lower.distT(2,:);

for i = 1:length(distL)-1
    i1          = (upper.distT(1,:) >= distL(1,i));
    i2          = (upper.distT(1,:) <= distL(1,i+1));
    index       = i1 == i2;
    uT(index)   = max(upper.distT(2,index));
    lT(index)   = max(lower.distT(2,index));
end

wing.t_s = upper.t_s;
wing.b = upper.b;
wing.h = upper.h;

t=[];
t(1:4)= 7.38e-3
t(5:7)= 6.76e-3
t(8:11)= 6.35e-3
t(12:14)=  5.68e-3
t(15:18)= 5.23e-3
t(19:21)= 4.65e-3
t(22:23)= 3.85e-3
t(24:25)= 2.21e-3

eff_thicknesses = zeros(length(t), 1);

for j = 1:length(t)
    wing.t = t(j);
    eff_thicknesses(j) = smearedThickness(wing);
end


save("effectiveThicknesses.mat", "eff_thicknesses", "rib")


