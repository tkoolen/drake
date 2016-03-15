function doubleIntegratorCapturability(n)
if nargin < 1
  n = 0;
end
model = DoubleIntegrator(1);

T = 1;

R_diag = [1.6 1.6];

options.degree = 4;
options.do_backoff = false;
options.backoff_ratio = 1.02;
options.scale = 1;
options.scale_input = 1;
options.control_design = true;
options.korda_control_design = true;
options.beta = .1;
options.infinite_time = false;

% R_diag = 2 * ones(1, model.num_states);

goal_radius = 0.1;
target = @(x) goal_radius^2 - x'*x;
% target = [];

[Vsol,Wsol,u_sol] = nStepCapturabilitySOS(model, T, R_diag, target, n, options);

options.degree = 8;
options.beta = .01;
[V_inner, W_inner] = finiteTimeInnerApproximation(model,T,u_sol,R_diag,target,options);

end
