function L = VTAero(VT, L_total)
% Takes VT (1 x n discretisation of y coordinates) and lift (OEI)
% Outputs L (1 x n array of sectional lift values) for load factor of 1
    

    index       = length(VT);
    s           = VT(index);

    L_0         = L_total / (s*pi/4)
    
    L           = L_0 * sqrt(1-(VT/s).^2).*(VT(2) - VT(1));
end