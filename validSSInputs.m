function valid = validSSInputs(geometry, material)
% calculates whether inputs are valid
% returns 1 if inputs are valid
% returns 0 if inputs are invalid

% initialise valid
valid = 1;

%% Compressive Buckling Coefficient
% K graph inputs
h_b     = geometry.h/geometry.b;
t_s_t   = geometry.t_s/geometry.t;

% check h_b is on the graph
if (h_b < 0.15) || (h_b > 1)
    valid = -1
elseif geometry.t_s < 0.001
    valid = false;
elseif ((t_s_t < 0.5) || (t_s_t > 2))
    valid = -2
else
    %% K
    K = compBuckCoef(h_b, t_s_t)
    %% Buckling Parameters
    I           = inertiaZ(geometry);
    trueN       = realN(geometry,0);
    
    %% Skin Buckling
    if plateBuckle(K, material.E, geometry) < trueN/geometry.t_sm
        valid = -3
    end
    
    %% Stringer Buckling
    if eulerBuckle(material.E, I, geometry) < trueN/geometry.t_sm
        valid = -4
    end
    
    %% Compressive Yield
    if material.compYield < trueN/geometry.t_sm
        valid = -5
    end
    
%     if valid > 0 
%         disp([plateBuckle(K, material.E, geometry), eulerBuckle(material.E, I, geometry), trueN/geometry.t_sm]);
%         valid = true;
%     end
end

if valid < 0
    valid = false;
elseif valid > 0 
    valid = true;
end