function tMod = modSkin(geometry, material, distT, sigmaR, distL)
% Takes discretised wing span and outputs skin thickness to resist combined
% loading

% Local shear
K_s = 4.85;
q = 1.56e6; %%%%%%%% DUMMY VALUE
tMod = zeros(size(distT));

for i = 1:length(distT)
    N = 3.75*(aeroMom(geometry, distT(i)) + Remerz(distT(i)));
    q = abs(Spar_z(geometry, material, distT(i)));
    
    %
    [~,idx] = min(abs(distT(i)-distL(1,:)));
    L = distL(2,idx);
    
    syms tShear
    eqn1 = N*geometry.b^2/(sigmaR*(tShear+geometry.A_s/geometry.b)*3.62*...
        material.E*tShear^2) + (q*geometry.b^2/(tShear^3*K_s*material.E))^2 == 0.99;
    tComb1 = double(vpasolve(eqn1,tShear));
    
    tComb1(real(tComb1)<0) = -inf;
    tComb1(imag(tComb1)~=0) = -inf;
    tMod1 = max(tComb1);
    
    % Global shear
    eqn2 = N*L^2/(sigmaR*(tShear+geometry.A_s/geometry.b)*3.62*...
        material.E*(tShear+geometry.A_s/geometry.b)^2) + ...
        (q*L^2/((tShear+geometry.A_s/geometry.b)^3*K_s*material.E))^2 == 0.99;
    tComb2 = double(vpasolve(eqn2,tShear));
    
    tComb2(real(tComb2)<0) = -inf;
    tComb2(imag(tComb2)~=0) = -inf;
    tMod2 = max(tComb2);
    
    tMod(i) = max([tMod1, tMod2]);
    
    if tMod(i) < 0.001
        tMod(i) = 0.001;
    end

end