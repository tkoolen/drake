function lipm3DNstepCapturability(n)

model_params.T = 0.3;  % Step-time
model_params.g = 10; % gravity acceleration
model_params.z_nom = 1; % nominal center of mass height
model_params.step_max = .7; % max step distance

options.filename_suffix = 'LIPM';
if n > 0
  options.time_scaling = model_params.T;
else
  options.time_scaling = 2;
end
options.plotfun = @(n, Vsol, Wsol, h_X, R_diag, t, x) lipm3DPlotFun(n, Vsol, Wsol, h_X, R_diag, t, x, model_params);

dynamics = @(t, x) lipmDynamics(t, x, model_params);
reset = @lipmResetMap;
reset_input_limits = @(u) lipmStepLimits(u, model_params);
nStepCapturabilitySOS(dynamics, reset, reset_input_limits, n, options)

end

function xdot = lipmDynamics(t, x, model_params)
g = model_params.g;
z_nom = model_params.z_nom;

q = x(1 : 2);
v = x(3 : 4);
xdot = [v; q*g/z_nom];
end

function xp = lipmResetMap(t, x, u)
% control input changes q only
% qp = qm - u
q = x(1 : 2);
v = x(3 : 4);
xp = [q - u; v];
end

function ret = lipmStepLimits(u, model_params)
step_max = model_params.step_max;
ret = step_max^2 - u'*u;
end
