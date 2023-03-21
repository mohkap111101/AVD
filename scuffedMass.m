function mass = scuffedMass(geometry, material, y)
% calculates mass of cover cross section

width   = geometry.wingboxWidth * geometry.C_r*(1-(1-geometry.taper)*(geometry.TEfus + y)/geometry.span);
A       = width * geometry.t;
n       = floor(width/geometry.b);
mass    = (n*geometry.A_s + A) * material.density;