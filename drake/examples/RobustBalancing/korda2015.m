function out = korda2015(model, options)

% dependencies
checkDependency('spotless');
checkDependency('mosek');

% options
R_diag = options.R_diag;
M = options.M;
ubar = options.ubar;
beta = options.beta;
d = options.degree;

% problem setup
x = msspoly('x', model.num_states);
g_X = 1 - x' * diag(1./(R_diag.^2)) * x;
l_x = options.l_x(x);
l_u = options.l_u(x);
if options.split_inputs
  l_u = [l_u; l_u];
end

[umin, umax] = model.simpleInputLimits(); % TODO: check that input limits are simple
if options.split_inputs
  u_0_to_ubar = msspoly('u', 2 * model.num_inputs);
  
  u_0_to_ubar_pos = u_0_to_ubar(1 : model.num_inputs);
  u_0_to_ubar_neg = u_0_to_ubar(model.num_inputs + 1 : end);
  
  u_pos = u_0_to_ubar_pos / ubar .* umax;
  u_neg = u_0_to_ubar_neg / ubar .* -umin;
  
  u = u_pos - u_neg;
else
  u_0_to_ubar = msspoly('u', model.num_inputs);
  u = umin + u_0_to_ubar / ubar .* (umax - umin);
end

xdot = model.dynamics(0, x, u);
G = diff(xdot, u_0_to_ubar); % TODO: check control affine
f = xdot - G * u_0_to_ubar;

prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog, rho] = prog.newFreePoly(monomials(x, 0 : d));
[prog, rho0] = prog.newFreePoly(monomials(x, 0 : d));
[prog, rhoT] = prog.newFreePoly(monomials(x, 0 : d));

m = size(G, 2);
sigma = zeros(m, 1) * msspoly(0);
for i = 1 : m
  [prog, sigma(i)] = prog.newFreePoly(monomials(x, 0 : d));
end

cost = spotlessIntegral(prog, l_x * rho + l_u' * sigma + M * rhoT, x, R_diag, [], []);
prog = prog.withPolyEqs(rhoT - rho0 + beta * rho + sum(sum(diff(rho * f, x))) + sum(sum(diff(G * sigma, x))));

gbar = prod(g_X);
[prog, rho_sos] = spotless_add_eq_sprocedure(prog, -rho, -gbar, x, d - even_degree(gbar, x));
for i = 1 : length(g_X)
  [prog, rho_sos] = spotless_add_eq_sprocedure(prog, rho_sos, -g_X(i), x, d - even_degree(g_X(i), x));
end
prog = prog.withSOS(rho_sos);

prog = prog.withSOS(rho0 - 1);

for i = 1 : m
  [prog, sigma_limit_sos] = spotless_add_sprocedure(prog, ubar * rho - sigma(i), gbar, x, d - even_degree(gbar, x));
  prog = prog.withSOS(sigma_limit_sos);
end

prog = prog.withSOS(rhoT);

for i = 1 : m
  [prog, sigma_sos] = spotless_add_sprocedure(prog, sigma(i), gbar, x, d - even_degree(gbar, x));
  prog = prog.withSOS(sigma_sos);
end

% solve
spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(cost, solver, spot_options);

out.x = x;
out.g_X = g_X;
out.sigma = sol.eval(sigma);
out.rho = sol.eval(rho);
out.f = f;
out.G = G;
out.fbar = @(x_val) closedLoopDynamics(f, G, out.sigma, out.rho, x, x_val);

end

function xdot = closedLoopDynamics(f, G, sigma, rho, x, x_val)
f = dmsubs(f, x, x_val);
G = reshape(dmsubs(G(:), x, x_val), size(G));
sigma = dmsubs(sigma, x, x_val);
rho = dmsubs(rho, x, x_val);
xdot = f + G * (sigma ./ rho);
end
