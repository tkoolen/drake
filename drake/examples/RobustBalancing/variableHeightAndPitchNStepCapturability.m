function variableHeightAndPitchNStepCapturability(n)

g = 10;
z_nom = 1;
R_diag = [2, .5, 1, 2, 2, 2];
f_min = .5;
f_max = 1.5;
inertia_ratio = .3^2/2; 
step_max = .7;
step_time = .3;

model = VariableHeightandPitch2D(g, z_nom, step_max, step_time, inertia_ratio, f_max, f_min);

options.degree = 4;
options.do_backoff = false;
options.backoff_ratio = 1.05;
options.free_final_time = false;
options.scale = 1/4;
options.scale_input = 2;

target = [];
% goal_radius = 0.01;
% target = @(x) goal_radius^2 - x'*x;

if n > 0
  T = step_time;
else
  T = 1;
end

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
