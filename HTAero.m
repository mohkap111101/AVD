function L = HTAero(HT, L_total)
% Takes HT (1 x n discretisation of y coordinates) and lift (maximum takeoff weight)
% Outputs L (1 x n array of sectional lift values) for load factor of 1
    
    fuselage    = 0.854;
    HT          = HT + fuselage;
    index       = length(HT);
    s           = HT(index);

    L_0         = L_total / (s*0.696);
    
    L           = L_0 * sqrt(1-(HT/s).^2).* (HT(2) - HT(1));
end