function mass = wingMass(geometry, material)
% calculates mass of cover cross section

geometry.t_sm = smearedThickness(geometry);

width   = geometry.wingboxWidth * geometry.C_r*(1-(1-geometry.taper)*geometry.TEfus/geometry.span);
A_s     = 1.6 * geometry.h * geometry.t_s * floor(width/geometry.b);
A       = width * geometry.t;
mass    = (A_s + A) * material.density;

if ~validSSInputs(geometry, material)
    mass = mass * 1e10;
end

end