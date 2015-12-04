function lipm3DNstepCapturability(n)

g = 10;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = 0.1; % set to 0 to get point foot model with no continuous inputs

model = LIPM3D(g, z_nom, step_max, step_time, cop_max);

if n > 0
  T = step_time;
else
  T = 2;
end
options.degree = 6;

R_diag = 2 * ones(1, model.num_states);

% radius of ball around the origin used as goal for 0-step capturability
goal_radius = 0.01;
target = @(x) goal_radius^2 - x'*x;

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
