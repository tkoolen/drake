function firstOrderSwingVariableHeightNStepCapturability(n)

g = 10;
u_max = 1;
z_nom = 1;
R_diag = [2, 1, 2, 2, 2];
f_min = .5;
f_max = 1.5;

model = FirstOrderSwingVariableHeight(g, z_nom, u_max, f_max, f_min);

options.degree = 4;
options.do_backoff = false;
options.backoff_ratio = 1.05;
options.free_final_time = true;
options.scale = 1/2;
options.scale_input = 2;
options.control_design = true;

target = [];
% goal_radius = 0.01;
% target = @(x) goal_radius^2 - x'*x;

if n == 0
  T = 1;
else
  T = .3;
end

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
