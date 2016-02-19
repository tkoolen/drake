function lipm2DNstepCapturability(n)
if nargin < 1
  n = 0;
end
g = 1;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = .1; % set to 0 to get point foot model with no continuous inputs

model = LIPM2D(g, z_nom, step_max, step_time, cop_max);
R_diag = [1 1];
if n > 0
  T = step_time;
else
  T = 5;
end
options.degree = 8;
options.scale = 1;
options.control_design = true;
% options.free_final_time = true;

% radius of ball around the origin used as goal for 0-step capturability
% goal_radius = 0.01;
% target = @(x) goal_radius^2 - x'*x;
target = [];

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
