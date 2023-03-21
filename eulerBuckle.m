function sigma = eulerBuckle(E, I, geometry)
% calculates whether Euler buckling occurs
% returns 1 if stringers buckle
% returns 0 if it is safe

A_s     = 1.6 * geometry.t_s * geometry.h;

P_crit  = pi^2 * E * I / geometry.L^2;

sigma   = P_crit/A_s;

end