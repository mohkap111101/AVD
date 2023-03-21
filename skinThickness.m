function t = skinThickness(geometry)
% Takes lifting surface parameters
% Calculates ideal skin thickness


% For Farrar efficiency = 0.95:
%       - T = t + A_s/b
%       - A_s/bt ~ 1.5
%       - t_s/t ~ 1.05
%       - sigma_crit/sigma_0 ~ 
t   = 2/5 * geometry.T;
A_s = 3/5 * geometry.T/geometry.pitch_s;