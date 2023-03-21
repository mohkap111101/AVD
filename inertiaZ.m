function I = inertiaZ(geometry)
% Calculates lesser inertia (I_xx or I_yy) for given stringer dimensions

Ixx = 7/30 * geometry.h^3 * geometry.t_s;
Iyy = 387/2000 * geometry.h^3 * geometry.t_s + geometry.t_s^3 * geometry.h / 12;
Ixy = 9/200 * geometry.h^3 * geometry.t_s;

theta = 0.5*atan(-2*Ixy/(Ixx-Iyy));
Ixmod = (Ixx + Iyy)/2 + (Ixx - Iyy)/2*cos(2*theta) - Ixy*sin(2*theta);
Iymod = (Ixx + Iyy)/2 - (Ixx - Iyy)/2*cos(2*theta) + Ixy*sin(2*theta);

if Ixx < Iyy
    I = Ixx;
else
    I = Iyy;
end
I = Iymod;
end