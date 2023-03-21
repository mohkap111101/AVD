clear
clc

b =  5.584776594*2;
z = 0:0.1:floor(b*10/2)/10;
L_req = 28.685 * 10^3;
L_0 = (8*L_req/(b*pi));
dL = L_0.*sqrt(1-(z./(b/2)).^2);
L = zeros(1,length(z));


for i = 1:length(z)-1
    
    L(i) = ((dL(i)+dL(i+1))/2)*(z(i+1)-z(i));
    
end

SF = zeros(size(L));

for i = [1:length(L)]

    SF(i) = sum(L([i:length(L)]));

end

dM = zeros(size(SF));

for i = [1:length(z)-1]
    
    dM(i) = ((SF(i)+SF(i+1))/2)*(z(i+1)-z(i));
    
end

BM = zeros(size(dM));

for i = [1:length(dM)]

    BM(i) = sum(dM([i:length(dM)]));

end




figure(1)
plot(z,dL, "-ro", "MarkerFaceColor", "r", "MarkerSize", 3)
title("Lift Distribution on Vertical tail")
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Lifting Load, N", "Interpreter","latex", "FontSize", 16)
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure(2)
plot(z,SF, "-ro", "MarkerFaceColor", "r", "MarkerSize", 3);
title("Shear Flow Distribution on Vertical tail")
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Shear Flow, N", "Interpreter","latex", "FontSize", 16)
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")

figure(3)
plot(z,BM, "-ro", "MarkerFaceColor", "r", "MarkerSize", 3);
title("Bending Moment Distribution on Vertical tail")
xlabel("Station along wing, m", "Interpreter","latex", "FontSize", 16)
ylabel("Bending Moment, Nm", "Interpreter","latex", "FontSize", 16)
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")





