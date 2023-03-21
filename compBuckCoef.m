function K = compBuckCoef(h_b, ts_t)
% Takes ts/t and h/b as inputs
% ts/t must be equal to one of [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5,
% 2.0]
% h/b should be between 0.15 and 1
% returns compressive buckling coefficient for z stringers

% Load coeffs
load('catchpoleCoeffs.mat');

ts_t_array = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2];

[idx, idx_1] = GetArrayIndicies(ts_t, ts_t_array);

% Get K for idx
if ts_t_array(idx) == 0.5
    K_idx = polyval(C05, h_b);
elseif ts_t_array(idx) == 0.6
    K_idx = polyval(C06, h_b);
elseif ts_t_array(idx) == 0.7
    K_idx = polyval(C07, h_b);
elseif ts_t_array(idx) == 0.8
    K_idx = polyval(C08, h_b);
elseif ts_t_array(idx) == 0.9
    K_idx = polyval(C09, h_b);
elseif ts_t_array(idx) == 1.0
    K_idx = polyval(C10, h_b);
elseif ts_t_array(idx) == 1.25
    K_idx = polyval(C125, h_b);
elseif ts_t_array(idx) == 1.5
    K_idx = polyval(C150, h_b);
elseif ts_t_array(idx) == 2.0
    K_idx = polyval(C20, h_b);
else
    K_idx = NaN;
end

% Get K for idx_1
if ts_t_array(idx_1) == 0.5
    K_idx_1 = polyval(C05, h_b);
elseif ts_t_array(idx_1) == 0.6
    K_idx_1 = polyval(C06, h_b);
elseif ts_t_array(idx_1) == 0.7
    K_idx_1 = polyval(C07, h_b);
elseif ts_t_array(idx_1) == 0.8
    K_idx_1 = polyval(C08, h_b);
elseif ts_t_array(idx_1) == 0.9
    K_idx_1 = polyval(C09, h_b);
elseif ts_t_array(idx_1) == 1.0
    K_idx_1 = polyval(C10, h_b);
elseif ts_t_array(idx_1) == 1.25
    K_idx_1 = polyval(C125, h_b);
elseif ts_t_array(idx_1) == 1.5
    K_idx_1 = polyval(C150, h_b);
elseif ts_t_array(idx_1) == 2.0
    K_idx_1 = polyval(C20, h_b);
else
    K_idx_1 = NaN;
end

% Linearly interpolate between the two data points of (ts_t at idx, K_idx)
% and (ts_t at idx_1, K_idx_1) to get K
K = interp1([ts_t_array(idx) ts_t_array(idx_1)], [K_idx K_idx_1], ts_t);
