function lipm3DNstepCapturability(n)
if nargin < 1
  n = 0;
end
g = 10;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = 0.1; % set to 0 to get point foot model with no continuous inputs

model = LIPM3D(g, z_nom, step_max, step_time, cop_max);
R_diag = [.5 .5 2 2];
if n > 0
  T = step_time;
else
  T = 1;
end
options.degree = 4;
options.scale = 1;
options.control_design = true;
% options.free_final_time = true;

% radius of ball around the origin used as goal for 0-step capturability
% goal_radius = 0.01;
% target = @(x) goal_radius^2 - x'*x;
target = [];

nStepCapturabilitySOS(model, T, R_diag, target, n, options)

if n == 0 && options.control_design
  sol = load(['V0_' class(model) '.mat']);
  t = sol.t;
  x = sol.x;

  test_manual_controller = true;
  if test_manual_controller
    q = x(1 : 2);
    v = x(3 : 4);
    omega0 = sqrt(g / z_nom);
    ric = q + v / omega0;
    k = 3; % >= 1 works
    u = k * ric / model.cop_max;
  else
    u = clean(sol.u_sol);
  end

  B = barrierFunctionForClosedLoopSystem(model, R_diag, t, x, u, struct('B_degree', 4));
  figure(5);
  contourSpotless(B, x(1), x(3), [-1 1], [-1 1], [x(2); x(4); t], [0; 0; 0], 0);
  xlabel('q_1'); ylabel('v_1'); legend('B(x)');
end

end
