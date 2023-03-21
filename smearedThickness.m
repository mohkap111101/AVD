function t_sm = smearedThickness(geometry)
% calculates smeared thickness for input parameters

A_s     = 1.6 * geometry.t_s * geometry.h;
t_sm    = geometry.t + A_s/geometry.b;

end