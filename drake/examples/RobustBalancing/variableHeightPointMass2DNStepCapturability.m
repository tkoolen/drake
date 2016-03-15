function variableHeightPointMass2DNStepCapturability(n)

if nargin < 1
  n = 0;
end

g = 3;
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
  T = 2;
end
options.degree = 6;
options.do_backoff = false;
options.backoff_ratio = 1.02;
options.scale = 1;
options.scale_input = 1;
options.control_design = true;
options.korda_control_design = true;
options.beta = .01;
options.infinite_time = false;

% R_diag = 2 * ones(1, model.num_states);

% goal_radius = 0.1;
% target = @(x) goal_radius^2 - x'*x;
target = [];

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
