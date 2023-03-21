function [SF, BM] = CalculateSFandBM(loads, y)
% Function to determine the shear force and bending moment distributions
% along the span on the wing given a particular load distribution

ds = y(2) - y(1);

SF = zeros(1, length(y));
d_BM = zeros(1, length(y));
BM = zeros(1, length(y));

for i=1:length(y)
    SF(i) = sum(loads(i:end));
end

% for i=1:length(y)-1
%     d_BM(i) = (SF(i+1) + SF(i))*ds/2;
% end
% 
% d_BM(end) = 0;

for i=1:length(y)
    BM(i)=0;
    for j=i:length(y)
        BM(i)=BM(i)+loads(j)*(y(j)-y(i));
    end
end