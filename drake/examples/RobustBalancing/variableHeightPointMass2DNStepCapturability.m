function variableHeightPointMass2DNStepCapturability(n)

g = 2;
step_max = .7;
step_time = 0.3;
u_max = 1.5 * g;
z_nom = 1;

model = VariableHeightPointMass2D(g, z_nom, step_max, step_time, u_max);

if n > 0
  T = step_time;
else
  T = 2;
end
% options.plotfun = @(n, Vsol, Wsol, h_X, R_diag, t, x) lipm3DPlotFun(n, Vsol, Wsol, h_X, R_diag, t, x, model);
options.degree = 6;

% R_diag = 2 * ones(1, model.num_states);
R_diag = [2, 0.5, 2, 2];

goal_radius = 0.01;
target = @(x) goal_radius^2 - x'*x;

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
