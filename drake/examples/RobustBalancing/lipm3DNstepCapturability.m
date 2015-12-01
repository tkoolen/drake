function lipm3DNstepCapturability(n)

g = 10; % gravitational acceleration
z_nom = 1; % nominal center of mass height
step_max = .7; % max step distance
step_time = 0.3;  % step time

model = LIPM3D(g, z_nom, step_max, step_time);

if n > 0
  T = step_time;
else
  T = 2;
end
options.plotfun = @(n, Vsol, Wsol, h_X, R_diag, t, x) lipm3DPlotFun(n, Vsol, Wsol, h_X, R_diag, t, x, model);
options.degree = 6;

R_diag = 2 * ones(1, model.num_states);

% radius of ball around the origin used as goal for 0-step capturability
goal_radius = 0.01;
target = @(x) goal_radius^2 - x'*x;

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
