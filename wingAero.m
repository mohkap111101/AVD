function L = wingAero(wing, MTOW)
% Takes wing (1 x n discretisation of y coordinates) and MTOW (maximum takeoff weight)
% Outputs L (1 x n array of sectional lift values) for load factor of 1
    
    fuselage    = 3.34/2;
    wing        = wing + fuselage;
    weight      = MTOW.*9.81; % MTOW (N)
    index       = length(wing);
    s           = wing(index);

    L_0         = (weight/2) / (s*0.672);
    
    L           = L_0 * sqrt(1-(wing/s).^2) .* (wing(2) - wing(1));
end