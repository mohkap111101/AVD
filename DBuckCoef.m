function K = DBuckCoef(x, a_b)

% Load coeffs
load('DSectionFn.mat');

a_b_array = [1 , 1.5 , 2 , 3];

[idx, idx_1] = GetArrayIndicies(a_b, a_b_array);

% Get K for idx
if a_b_array(idx) == 1
    K_idx = polyval(One, x);
elseif a_b_array(idx) == 1.5
    K_idx = polyval(OnePointFive, x);
elseif a_b_array(idx) == 2
    K_idx = polyval(Two, x);
elseif a_b_array(idx) == 3
    K_idx = polyval(Three, x);
else
    K_idx = NaN;
end

% Get K for idx_1
if a_b_array(idx_1) == 1
    K_idx_1 = polyval(One, x);
elseif a_b_array(idx_1) == 1.5
    K_idx_1 = polyval(OnePointFive, x);
elseif a_b_array(idx_1) == 2
    K_idx_1 = polyval(Two, x);
elseif a_b_array(idx_1) == 3
    K_idx_1 = polyval(Three, x);
else
    K_idx_1 = NaN;
end

% Linearly interpolate between the two data points of (ts_t at idx, K_idx)
% and (ts_t at idx_1, K_idx_1) to get K
K = interp1([a_b_array(idx) a_b_array(idx_1)], [K_idx K_idx_1], a_b);