function firstOrderSwingZMPNStepCapturability(n)

g = 10;
u_max = 1;
T = .5;
z_nom = 1;
R_diag = [2, 2, 2];

scale = 1;
g = g/scale;
z_nom = z_nom/scale;
R_diag = R_diag/scale;

model = FirstOrderSwingZMP(g, z_nom, u_max);

options.degree = 4;
options.do_backoff = false;
options.backoff_ratio = 1.05;
options.free_final_time = true;
options.control_design = true;

target = [];
% goal_radius = 0.01;
% target = @(x) goal_radius^2 - x'*x;

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
