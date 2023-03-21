function K = DBuckCoefb(x, b_a)

% Load coeffs
load('DSectionFnb.mat');

b_a_array = [1 , 1.5 , 5];
[idx, idx_1] = GetArrayIndicies(b_a, b_a_array);

% Get K for idx
if b_a_array(idx) == 1
    K_idx = polyval(bOne, x);
elseif b_a_array(idx) == 1.5
    K_idx = polyval(bOnePointFive, x);
elseif b_a_array(idx) == 5
    K_idx = polyval(bFive, x);
else
    K_idx = NaN;
end

% Get K for idx_1
if b_a_array(idx_1) == 1
    K_idx_1 = polyval(bOne, x);
elseif b_a_array(idx_1) == 1.5
    K_idx_1 = polyval(bOnePointFive, x);
elseif b_a_array(idx_1) == 5
    K_idx_1 = polyval(bFive, x);
else
    K_idx_1 = NaN;
end

% Linearly interpolate between the two data points of (ts_t at idx, K_idx)
% and (ts_t at idx_1, K_idx_1) to get K
K = interp1([b_a_array(idx) b_a_array(idx_1)], [K_idx K_idx_1], b_a);