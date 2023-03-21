function [t1, b, A] = GetShomeGoodShparsh(t2, h, I_required)

% Setting up the optimisation problem
Area = @(t1, b) 2*b*t1 + (h - 2*t1)*t2;
Mass = @(tf, hf, bf) Area(tf, hf, bf) * 2*pi*Fuselage.cabin_radius * Fuselage.frame_density*(cabin_length./Lf);
find_tf = @(hf, bf, Lf) I_required(Lf)./(hf.^3/12 + 2*bf.*(hf./2).^2);
obj = @(tf, hf, bf) Area(teq(tf, hf, bf), hf, bf);
% Defining the non linear I constraint

c = @(x)[] ;%sigma_a - Fuselage.stringer_E*Is*pi^2/(x(4)^2 * As);
ceq = @(x) (1/6 * x(2) * x(1)^3 + 1/2 * x(2)* x(1) * (h - x(1))^2 + 1/12 * t2 * (h - 2*x(1))^3) - I_required;

Iconstraint = @(x)deal(c(x), ceq(x));

t1_min = 0.001; % 1mm minimum thickness
t1_max = 0.10; % 10cm max thickness

b_min = 2e-2 + t2;
b_max = 10e-2;


X0 = [3, 3].*(10^-2);

options = optimoptions('fmincon','Display','iter', "Algorithm","interior-point",...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg");
[minArea_dims, A] = fmincon(@(x) Area(x(1), x(2)), X0, [], [], [], [], [t1_min, b_min], [t1_max, b_max], Iconstraint, options);

t1 = minArea_dims(1);
b = minArea_dims(2);



end