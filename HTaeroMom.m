function BM = HTaeroMom(geometry, z)
% Takes wing (1 x n discretisation of y coordinates) and MTOW (maximum takeoff weight)
% Outputs L (1 x n array of sectional lift values) for load factor of 1
    
    weight      = 33000; % MTOW (N)
    disc        = linspace(0,geometry.span/2,1000);

    L_0         = weight / (2*geometry.span*0.672);
    L           = L_0 * sqrt(1-(disc/(geometry.span/2)).^2);
    
    index       = (disc-geometry.TEfus)>z;
    if all(index(1:end-1)==0)
        BM = 0;
    else
        BM          = trapz(disc(index),L(index).*(disc(index)-z));
    end
end