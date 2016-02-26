function out = korda2014RegionOfAttractionAndController(model, target, options)

% dependencies
checkDependency('spotless');
checkDependency('mosek');

% options
v_degree = options.v_degree;
betas = options.betas;
betas_outer = betas(options.beta_outer_ind);
R_diag = options.R_diag;

% shared problem setup
ubar = 1;
x = msspoly('x', model.num_states);
g_X = 1 - x' * diag(1./(R_diag.^2)) * x;

u_0_to_ubar = msspoly('u', model.num_inputs);
[umin, umax] = model.simpleInputLimits(); % TODO: check that input limits are simple
u = umin + u_0_to_ubar / ubar * (umax - umin);
xdot = model.dynamics(0, x, u);
G = diff(xdot, u_0_to_ubar); % TODO: check control affine
f = xdot - G * u_0_to_ubar;

out = struct;
out.x = x;
out.g_X = g_X;
[outer_sol, out.v_outer, mu_ind, sigma_inds] = regionOfAttractionOuterApprox(f, G, ubar, g_X, R_diag, target, v_degree, betas_outer);
[out.u_0_to_ubar, out.u_hat_0_to_ubar] = extractController(outer_sol, ubar, g_X, mu_ind, sigma_inds);
out.u = (out.u_0_to_ubar - umin) / (umax - umin) * ubar;
out.u_hat = (out.u_hat_0_to_ubar - umin) / (umax - umin) * ubar;
out.fbar = f + G * out.u_0_to_ubar; % closed loop dynamics
out.fbar_hat = f + G * out.u_hat; % closed loop dynamics for controller without input limits
out.vs_inner = regionOfAttractionInnerApprox(out.fbar, g_X, R_diag, target, betas, options.w_v_degree);

end

function [sol, v_sol, mu_ind, sigma_inds] = regionOfAttractionOuterApprox(f, G, ubar, g_X, R_diag, target, v_degree, beta)
% sos program setup
prog = spotsosprog;
x = decomp(g_X);
prog = prog.withIndeterminate(x);
[prog, v] = prog.newFreePoly(monomials(x, 0 : v_degree));

dynamics_sos = beta * v - diff(v, x) * f;
p_degree = even_degree(dynamics_sos, x);
m = size(G, 2);
sigma_inds = zeros(m, 1);
for i = 1 : m
  [prog, p_i] = prog.newFreePoly(monomials(x, 0 : p_degree));
  
  p_i_affine_part_sos = p_i - diff(v, x) * G(:, i);
  [prog, p_i_affine_part_sos] = spotless_add_sprocedure(prog, p_i_affine_part_sos, g_X, x, p_degree - 2);
  [prog, sigma_inds(i)] = prog.withSOS(p_i_affine_part_sos);
  
  [prog, p_i_sos] = spotless_add_sprocedure(prog, p_i, g_X, x, p_degree - 2);
  prog = prog.withSOS(p_i_sos);  
  
  dynamics_sos = dynamics_sos - ubar * p_i;
end
[prog, dynamics_sos] = spotless_add_sprocedure(prog, dynamics_sos, g_X, x, even_degree(dynamics_sos, x) - 2);
[prog, mu_ind] = prog.withSOS(dynamics_sos);

v_target_sos = v - 1;
if isnumeric(target)
  prog = prog.withPos(subs(v_target_sos, x, target));
else
  [prog, v_target_sos] = spotless_add_sprocedure(prog, v_target_sos, g_X_target, x, v_degree - 2);
  prog = prog.withSOS(v_target_sos);
end

v_sos = v + 1;
[prog, v_sos] = spotless_add_sprocedure(prog, v_sos, g_X, x, v_degree - 2);
prog = prog.withSOS(v_sos);

cost = spotlessIntegral(prog, v, x, R_diag, [], []);

% solve
spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(cost, solver, spot_options);
v_sol = sol.eval(v);
xlabel('x_1')
ylabel('x_2')

end

function [u_sol, u_hat_sol] = extractController(sol, ubar, g_X, mu_ind, sigma_inds)

% extract controller that doesn't satisfy input limits
x = decomp(g_X);
k = even_degree(sol.prog.sosExpr(mu_ind), x) / 2;
y_mu = sol.prog.sosEqsDualVars{mu_ind}; % moments of mu
basis_y_mu = sol.prog.sosEqsBasis{mu_ind}; % TODO: check
basis_y_mu_grlex = grlex(basis_y_mu);
y_mu_to_grlex = match_monomials(basis_y_mu_grlex, basis_y_mu);
if ~all(double(basis_y_mu(y_mu_to_grlex) - basis_y_mu_grlex) == 0)
  error('logic mistake');
end
y_mu_grlex = y_mu(y_mu_to_grlex);
y_mu_grlex_sol = double(sol.dualEval(y_mu_grlex));
M = momentMatrix(y_mu_grlex_sol, basis_y_mu_grlex);
m = length(sigma_inds); % num_inputs
n = size(x, 1); % num_states
if nchoosek(n + k, n) ~= size(M, 1)
  error('logic mistake')
end

u_hat_sol = zeros(m, 1) * msspoly(0);
u_hat_coeffs = cell(m, 1);
for i = 1 : m
  y_sigma = sol.prog.sosEqsDualVars{sigma_inds(i)};
  basis_y_sigma = sol.prog.sosEqsBasis{sigma_inds(i)};
  basis_y_sigma_grlex = grlex(basis_y_sigma);
  y_sigma_to_grlex = match_monomials(basis_y_sigma_grlex, basis_y_sigma);
  if ~all(double(basis_y_sigma(y_sigma_to_grlex) - basis_y_sigma_grlex) == 0)
    error('logic mistake');
  end
  y_sigma_grlex = y_sigma(y_sigma_to_grlex);
  y_sigma_grlex_sol = double(sol.dualEval(y_sigma_grlex));
%   if ~all(double(grlex(sol.prog.gramMonomials{sigma_inds(i)}) - basis_y_sigma_grlex(1 : size(M, 1))) == 0)
%     error('logic mistake')
%   end
  % (14) in Korda, Henrion, Jones 2014
%   u_hat_coeffs{i} = M \ y_sigma_grlex_sol(1 : size(M, 1));
  
  M_sigma = momentMatrix(y_sigma_grlex_sol, basis_y_sigma_grlex);
  [U,S,V] = svd(M);
  s = diag(S);
  condTh = 1e6;
  keep = find((s(1)./s) < condTh,1,'last');
  v = U\M_sigma(:,1);
  v = diag(1./s(1:keep))*v(1:keep,:);
  u_hat_coeffs{i} = V(:,1:keep)*v;
  
  u_sol_basis = basis_y_sigma_grlex(1 : size(M, 1));
  u_hat_sol(i) = u_sol_basis' * u_hat_coeffs{i}; % controller violating input constraints
end

% make controller satisfy input limits
u_sol = zeros(m, 1) * msspoly(0);
for i = 1 : m
  % (15) in Korda, Henrion, Jones 2014
  u_degree = deg(u_sol_basis);
  M1_dim_and_M2_rows = nchoosek(n + u_degree, n);
  M2_cols = nchoosek(n + k, n);
  M1 = M(1 : M1_dim_and_M2_rows, 1 : M1_dim_and_M2_rows);
  M2 = M(1 : M1_dim_and_M2_rows, 1 : M2_cols);
  
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  [prog, u_i_coeffs] = prog.newFree(length(u_sol_basis), 1);
  u_i = u_i_coeffs' * u_sol_basis;
  %cost = u_i_coeffs' * M1 * u_i_coeffs - 2 * u_hat_coeffs{i}' * M2 * u_i_coeffs;
  %prog = prog.withPos(slack - cost); % doesn't work
  
  P = 2 * M1;
  lambda_min = min(eig(P));
  if lambda_min < 0
    if lambda_min > -1e-5
      P = P + 1.1 * abs(lambda_min) * eye(size(P));
    else
      error('P is not positive definite');
    end
  end
  q = -2 * (u_hat_coeffs{i}' * M2)';
  W = chol(P)';
  [prog, t] = prog.newFree(1, 1);
  cost = t + 2 * q' * u_i_coeffs;
  prog = prog.withPSD([eye(size(W)), W' * u_i_coeffs;
                       u_i_coeffs' * W, t]);
  
  [prog, u_lower_sos] = spotless_add_sprocedure(prog, u_i, g_X, x, even_degree(u_i, x) - 2);
  prog = prog.withSOS(u_lower_sos);
  
  [prog, u_upper_sos] = spotless_add_sprocedure(prog, ubar - u_i, g_X, x, even_degree(u_i, x) - 2);
  prog = prog.withSOS(u_upper_sos);
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(cost, solver, spot_options);
  
  u_sol(i) = sol.eval(u_i);
end

end

function vs_inner_sol = regionOfAttractionInnerApprox(fbar, g_X, R_diag, g_X_target, betas, w_v_degree)
prog = spotsosprog;
x = decomp(g_X);
prog = prog.withIndeterminate(x);
[prog, w] = prog.newFreePoly(monomials(x, 0 : w_v_degree));

if isnumeric(g_X_target)
%   g_X_target = -(x - g_X_target)' * (x - g_X_target);
  g_X_target = -g_X;
end

[prog, w_sos] = spotless_add_sprocedure(prog, w, -g_X_target, x, w_v_degree - 2);
prog = prog.withSOS(w_sos);

nbeta = length(betas);
vs = zeros(nbeta, 1) * msspoly(0);
w_v_sos = w - 1;
for i = 1 : nbeta
  [prog, v] = prog.newFreePoly(monomials(x, 0 : w_v_degree));
  vs(i) = v;
  
  vdot_sos = betas(i) * v - diff(v, x) * fbar;
  [prog, vdot_sos] = spotless_add_sprocedure(prog, vdot_sos, -g_X_target, x, even_degree(vdot_sos, x) - 2);
  prog = prog.withSOS(vdot_sos);
  
  w_v_sos = w_v_sos - v;
  
  [prog, v_sos] = spotless_add_eq_sprocedure(prog, v, g_X, x, w_v_degree - 2);
  prog = prog.withSOS(v_sos);
end
[prog, w_v_sos] = spotless_add_sprocedure(prog, w_v_sos, -g_X_target, x, w_v_degree - 2);
prog = prog.withSOS(w_v_sos);

cost = spotlessIntegral(prog, w, x, R_diag, [], []);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(cost, solver, spot_options);
vs_inner_sol = sol.eval(vs);
end
