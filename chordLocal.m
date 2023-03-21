function c = chordLocal(geometry, y)
% returns local geometry for given lifting surface geometry and spanwise
% location
c = geometry.C_r*(1-(1-geometry.taper)*(geometry.TEfus + y)/(geometry.span/2));