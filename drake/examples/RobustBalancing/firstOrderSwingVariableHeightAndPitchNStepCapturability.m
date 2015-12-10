function firstOrderSwingVariableHeightAndPitchNStepCapturability(n)

g = 10;
u_max = 1;
T = .5;
z_nom = 1;
R_diag = [2, 1, 1, 2, 2, 2, 1];
f_min = .5;
f_max = 1.5;
inertia_ratio = .3^2/2; 

model = FirstOrderSwingVariableHeightandPitch(g, z_nom, inertia_ratio, u_max, f_max, f_min);

options.degree = 4;
options.do_backoff = false;
options.backoff_ratio = 1.05;
options.free_final_time = true;

target = [];
% goal_radius = 0.01;
% target = @(x) goal_radius^2 - x'*x;

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
