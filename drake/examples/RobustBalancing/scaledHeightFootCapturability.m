function scaledHeightFootCapturability(n)
if nargin < 1
  n = 0;
end

g = 10;
step_max = .7;
step_time = 0.3;
z_nom = 1;
R_diag = [1, 1, 3, 3];

f_range = .4;
foot_radius = .5;

model = ScaledHeightFootModel(g, z_nom, step_max, step_time, f_range, foot_radius);

if n > 0
  T = step_time;
  options.free_final_time = true;
else
  T = 1;
%   options.free_final_time = true;
end
options.degree = 4;
options.do_backoff = false;
options.backoff_ratio = 1.02;
options.scale = 1./R_diag';
options.scale_input = 1;
options.control_design = true;
options.korda_control_design = true;
options.beta = 0;
options.infinite_time = false;

% R_diag = 2 * ones(1, model.num_states);

% goal_radius = 0.1;
% target = @(x) goal_radius^2 - x'*x;
target = [];
R_goal = [.1 .2 .1 .1];
A_goal = diag(1./R_goal.^2);
% target = @(x) 1 - x'*A_goal*x;

[Vsol,Wsol,u_sol] = nStepCapturabilitySOS(model, T, R_diag, target, n, options);

keyboard
[V_inner, W_inner] = finiteTimeInnerApproximation(model,T,u_sol,R_diag,target,options);
end
