function [SF,M] = SF_and_BM_daquing(loads,y)

ds = y(2) - y(1);

SF = zeros(1,length(y));
M = zeros(1,length(y));
dM = zeros(1,length(y));
BM = zeros(1,length(y));

for i = [1:length(y)]
    SF(i) = sum(loads(1:i)); 
end

for i = [2:length(y)]
    dM(i) = (1/2) * (SF(i) + SF(i-1)) * (y(i) - y(i-1));
end

for i = [1:length(y)]
    M(i) = sum(dM(1:i)); 
end


