function nonlinearDoubleIntegratorCapturability(n)
if nargin < 1
  n = 0;
end
model = NonlinearDoubleIntegrator();

T = 1;

R_diag = [1.2 1.2].^.5;

options.degree = 6;
options.do_backoff = false;
options.backoff_ratio = 1.02;
options.scale = 1;
options.scale_input = 1;
options.control_design = true;
options.korda_control_design = true;
options.beta = 1;
options.infinite_time = true;


% R_diag = 2 * ones(1, model.num_states);

goal_radius = 0.1;
target = @(x) goal_radius^2 - x'*x;
% target = [];

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

data = load('V0_NonlinearDoubleIntegrator');
options.beta = [10 1 .1 .01 .001];
keyboard
end
