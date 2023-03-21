function N = realN(geometry, y)
% Calculates N (compressive force per unit length in cover due to 
% bending moments

chord   = geometry.C_r*(1-(1-geometry.taper)*(y+geometry.TEfus)/(geometry.span/2));
M       = 3.75*(aeroMom(geometry, y) + Remerz(y));
N       = M/(geometry.wingboxWidth*geometry.maxHeight*chord);

end