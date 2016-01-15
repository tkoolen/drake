function variableHeightPointMass2DNStepCapturability(n)

g = 10;
step_max = .7;
step_time = 0.3;
z_nom = 1;
R_diag = [2, .5 2, 2];

f_max = 1.2;
f_min = .8;

% f_max = 2;
% f_min = 0;

model = VariableHeightPointMass2D(g, z_nom, step_max, step_time, f_max, f_min);

if n > 0
  T = step_time;
else
  T = 1;
end
options.degree = 6;
options.do_backoff = true;
options.backoff_ratio = 1.05;
options.scale = 1/4;
options.scale_input = 1;
options.control_design = true;

% R_diag = 2 * ones(1, model.num_states);

% goal_radius = 0.05;
% target = @(x) goal_radius^2 - x'*x;
target = [];

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
