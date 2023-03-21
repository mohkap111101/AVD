function sigma = plateBuckle(K, E, geometry)
% calculates critical buckling stress for the skin
% returns 1 if skin buckles
% returns 0 if it is safe

sigma   = K * E * (geometry.t/geometry.b)^2;

end