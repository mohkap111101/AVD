function L = HTfindL(geometry, material, I, N)
% Takes these variables and gives minimum rib spacing at station
safety      = 1.15;
Lstringer   = sqrt(pi^2 * material.E * I / (geometry.A_s*N/geometry.t_sm));
Lstringer   = 1000;
Lpanel      = sqrt(3.62*material.E*geometry.t_sm^2/(safety*N/geometry.t_sm));
L           = min([Lstringer, Lpanel]);