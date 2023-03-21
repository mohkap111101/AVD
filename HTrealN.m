function N = HTrealN(geometry, y)
% Calculates N (compressive force per unit length in cover due to 
% bending moments

chord   = geometry.C_r*(1-(1-geometry.taper)*(y+geometry.TEfus)/(geometry.span/2));
aeroM   = HTaeroMom(geometry, y);
[BM,~]  = HTweight(geometry, y);
BM      = 4.5 * BM;
M       = aeroM - BM;
N       = M/(geometry.wingboxWidth*geometry.maxHeight*chord);


end