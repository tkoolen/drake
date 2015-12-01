function lipm3DNstepCapturability(n)

g = 10; % gravitational acceleration
z_nom = 1; % nominal center of mass height
step_max = .7; % max step distance
T = 0.3;  % step time

model = LIPM3D(g, z_nom, step_max, T);

options.filename_suffix = 'LIPM';
if n > 0
  options.time_scaling = model.T;
else
  options.time_scaling = 2;
end
options.plotfun = @(n, Vsol, Wsol, h_X, R_diag, t, x) lipm3DPlotFun(n, Vsol, Wsol, h_X, R_diag, t, x, model);

nStepCapturabilitySOS(model, n, options)

end
