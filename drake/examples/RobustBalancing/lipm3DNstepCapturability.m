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
    k = 3;
    u = k * ric;
  else
    u = clean(sol.u_sol);
  end

  B = barrierFunctionForClosedLoopSystem(model, R_diag, t, x, u, struct('B_degree', 4));
  figure(4);
  contourSpotless(B, x(1), x(3), [-1 1], [-1 1], [x(2); x(4)], [0; 0], 0);
end

end

function B_sol = barrierFunctionForClosedLoopSystem(model, R_diag, t, x, u, options)

xdot_closed_loop = model.dynamics(t, x, u);

prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog, B] = prog.newFreePoly(monomials(x, 0 : options.B_degree));

% initial condition (and fix scaling of problem)
x0 = zeros(model.num_states, 1);
prog = prog.withEqs(subs(B, x, x0) + 1);

% Unsafe set constraint
A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x; % < 0 for failed states
[prog, B_failed_sos] = spotless_add_sprocedure(prog, B, -h_X, x, options.B_degree - deg(h_X));
prog = prog.withSOS(B_failed_sos); % - 1e-5);

% Barrier function derivative constraint
Bdot = diff(B, x) * xdot_closed_loop;
prog = prog.withSOS(-Bdot);

% Solve
spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = false;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(subs(B, x, x0), solver, spot_options);
B_sol = sol.eval(B);

end
