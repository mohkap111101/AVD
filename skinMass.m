function mass = skinMass(geometry, material, distT)
% Returns mass of skin based on geometry and thickness distribution

distT = [distT, [geometry.span/2 - geometry.TEfus; distT(2,end)]];
width = geometry.wingboxWidth*geometry.C_r*(1-(1-geometry.taper)*...
        (geometry.TEfus+distT(1,:))/(geometry.span/2));
mass = trapz(distT(1,:),distT(2,:).*width*material.density);

end