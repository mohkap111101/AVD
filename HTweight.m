function [BM, SF] = HTweight(geometry, y)
% Takes discretised horizontal tail coordinates (HT) and its total weight
% Outputs sectional weight at each section according to trapezoidal
% distribution of chord
weight_total = 360*9.81/2;
ht  = linspace(0,geometry.span/2, 1000);
w_r = 2 * weight_total/((1+geometry.taper)*geometry.span);
weight  = w_r * (1 - (1-geometry.taper) * (geometry.TEfus + ht)/(geometry.span/2));

index = (ht-geometry.TEfus) > y;

if all(index(1:end-1)==0)
    BM = 0;
    SF = 0;
else
    BM = trapz(ht(index),weight(index).*(ht(index)-geometry.TEfus));
    SF = trapz(ht(index),weight(index));
end


end