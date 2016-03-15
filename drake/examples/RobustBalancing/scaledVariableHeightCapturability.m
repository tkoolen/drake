function scaledVariableHeightCapturability(n)
if nargin < 1
  n = 0;
end

g = 10;
step_max = .7;
step_time = 0.3;
z_nom = 1;
R_diag = [1, 1, 1, 1];

f_range = .4;

model = ScaledVariableHeightModel(g, z_nom, step_max, step_time, f_range);

if n > 0
  T = step_time;
  options.free_final_time = true;
else
  T = 1;
  options.free_final_time = false;
end
options.degree = 6;
options.do_backoff = false;
options.backoff_ratio = 1.02;
options.scale = 1;%./R_diag';
options.scale_input = 1;
options.control_design = true;
options.korda_control_design = false;


% R_diag = 2 * ones(1, model.num_states);

% goal_radius = 0.1;
% target = @(x) goal_radius^2 - x'*x;
target = [];

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
