function [SF,M] = CalculateSFandBM_fuselage(loads,y)

ds = y(2) - y(1);

SF = zeros(1,length(y));
M = zeros(1,length(y));

for i = [1:length(y)]
    SF(i) = sum(loads(1:i)); 
end

for i = [2:length(y)]
    M(i) = M(i-1) + SF(i-1)*ds;
end

