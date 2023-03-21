function M_0 = wingMom(wing, M_0_total)
% Takes wing (1 x n discretisation of y coordinates) and M_0_total (total
% zero-lift pitching moment)
% Outputs M_0 (1 x n array of sectional zero-lift pitching moments)

    fuselage    = 3.1911/2;
    wing        = wing + fuselage;
    index       = length(wing);
    s           = wing(index);

    M_0_centre  = M_0_total/2 / (s*0.672);
    
    M_0         = M_0_centre * sqrt(1-(wing/s).^2);
end