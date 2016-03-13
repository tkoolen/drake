function orbitalEnergySym()
syms g real;
syms x z real;
syms xd zd real;
syms x0 z0 real;
syms xd0 zd0 real;
syms zf real;

q = [x; z];
v = [xd; zd];
q0 = [x0; z0];
v0 = [xd0; zd0];
n = 4;

% find trajectory
[fw, u_trajw, uw, w] = captureHeightTrajectory(g, q, v, q0, v0, zf, n);

% plug in g and zf
g_val = 9.8;
zf_val = 1;
fw = subs(fw, [g; zf], [g_val; zf_val]);
uw = subs(uw, [g; zf], [g_val; zf_val]);
u_trajw = subs(u_trajw, [g; zf], [g_val; zf_val]);

% plug in w
w_val = 0 * ones(size(w));
f = subs(fw, w, w_val);
u = subs(uw, w, w_val);
u_traj = subs(u_trajw, w, w_val);

% closed loop dynamics
vdot = [0; -g_val] + u * q;
xdot = simplify([v; vdot]);

% boundary function
[B, state, nu, du] = captureHeightControlSOS(xdot, u, q, v, g_val, zf_val);

z0_val = zf_val;
zd0_val = 0;

B_fig_num = 4;
figure(B_fig_num);
clf;
omega = sqrt(g_val / zf_val);
hold on
x_range = [-2, 2];
xd_range = [0, max(x_range) * omega];
% contourSpotless(subs(B, [state(2); state(4)], [z0_val; zd0_val]), state(1), state(3), x_range, xd_range, [], [], 0, {'b'});
contourSpotless(subs(nu, [state(2); state(4)], [z0_val; zd0_val]), state(1), state(3), x_range, xd_range, [], [], 0, {'b'});
% contourSpotless(subs(du, [state(2); state(4)], [z0_val; zd0_val]), state(1), state(3), x_range, xd_range, [], [], 0, {'r'});
xs = linspace(x_range(1), 0, 100);
icp_line = -omega * xs;
plot(xs, icp_line, 'r');
xlabel('q_x');
ylabel('v_x');

while true
  % evaluate resulting trajectory
  figure(B_fig_num);
  [x0_val, xd0_val, button] = ginput(1);
  if button == 27 % escape
    break;
  end
  
  % x0_val = -1;
  % omega0 = sqrt(g_val / z0_val);
  % xd0_val = -omega0 * x0_val * 0.95; %1.05;
  q0_val = [x0_val; z0_val];
  v0_val = [xd0_val; zd0_val];
  [xs, xds, feasible] = evaluateTrajectory(f, u_traj, u, xdot, q, v, q0, v0, q0_val, v0_val, g_val);
  
  figure(B_fig_num)
  hold on;
  if feasible
    color = 'g';
  else
    color = 'r';
  end
  plot(xs, xds, color);
  axis([x_range, xd_range]);
  hold off;
end

end

function [height_traj, u_traj, controller, w] = captureHeightTrajectory(g, q, v, q0, v0, zf, n)
x = q(1);
z = q(2);
xd = v(1);
zd = v(2);
x0 = q0(1);
z0 = q0(2);
xd0 = v0(1);
zd0 = v0(2);

a = sym('a', [n + 1, 1], 'real');

mons = sym.zeros(n + 1, 1);
for i = 0 : n
  mons(i + 1) = x^i;
end

f = a' * mons;
fp = diff(f, x);
fpp = diff(fp, x);

A = sym.zeros(4, n + 1);
A(1, :) = jacobian(subs(f, x, 0), a); % final height
A(2, :) = jacobian(subs(f, x, x0), a); % initial height
A(3, :) = jacobian(subs(fp, x, x0), a); % initial slope
for i = 0 : n
  A(4, i + 1) = 2 * g * (1 - i) / (i + 2) * x0^(i + 2); % orbital energy
end

b = [...
  zf;
  z0;
  zd0 / xd0;
  (xd0 * z0 - zd0 * x0)^2]; % same as xd0^2 * (z0 - zd0 / xd0 * x0)^2]; % note: xd0^2 is xd0 in paper, but that's wrong.

w = sym('w', [length(a) - size(A, 1), 1], 'real');

a_sol = A \ b + null(A) * w;

h = f - fp * x;
fint = int(f * x, x);
% orbital_energy = simplify(1/2 * xd^2 * h^2 + g * x^2 * f - 3 * g * fint);

height_traj = simplify(subs(f, a, a_sol));
u = (g + fpp * xd^2) / (z - fp * x);
u = simplify(subs(u, a, a_sol));
controller = simplify(subs(u, [x0; z0; xd0; zd0], [x; z; xd; zd]));

xd_squared = simplify(2 * (3 * g * fint - g * x^2 * f) / h^2);
u_traj = simplify((g + fpp * xd_squared) / (f - fp * x));
u_traj = simplify(subs(u_traj, a, a_sol));

% u_traj_den = simplify(h^2 * g + fpp * 2 * (3 * g * fint - g * x^2 * f));
% u_traj_num = simplify(h^3);

% uhsq_den = g * h^2 + 2 * g * fpp * (int(h * x, x) - subs(int(h * x, x), x, x0)) + fpp * xd0^2 * subs(h, x, x0)^2;

% uh_den

% bla = g * h / (2 * g * (int(h * x, x) - subs(int(h * x, x), x, x0)) + xd0^2 * subs(h^2, x, x0)) + fpp / h;

end

function xds = horizontalVelocities(f, g, xd0, xs)
x = symvar(f);
fp = diff(f, x);
h = f - fp * x;
fint = int(f * x, x);

xdsquared_numerator = 2 * (3 * g * fint - g * x^2 * f);
xdsquared_denominator = h^2;
xdsquared_numerator_val = polyval(sym2poly(xdsquared_numerator), xs);
xdsquared_denominator_val = polyval(sym2poly(xdsquared_denominator), xs);
xds = abs(sqrt(xdsquared_numerator_val ./ xdsquared_denominator_val)) * sign(xd0);
end

function f_legs = legForce(u, q, v, xs, zs, xds, zds)
fun = matlabFunction(u * norm(q), 'Vars', [q; v]);
f_legs = fun(xs, zs, xds, zds);
end

function ret = captureHeightTrajectorySOS(utraj_sym, x_sym, x0, u_max)
checkDependency('spotless');
checkDependency('mosek');

prog = spotsosprog;
[prog, x] = prog.newIndeterminate('x');

[n_sym, d_sym] = numden(utraj_sym);
n = sym2msspoly(x_sym, x, n_sym);
d = sym2msspoly(x_sym, x, d_sym);

% q0 = msspoly('y', length(q0_sym));
% v0 = msspoly('z', length(q0_sym));
% n = sym2msspoly([x_sym; q0_sym; v0_sym; w_sym], [x; q0; v0; w], n_sym);
% d = sym2msspoly([x_sym; q0_sym; v0_sym; w_sym], [x; q0; v0; w], d_sym);
% n = subs(n, [q0; v0], [q0_val; v0_val]);
% d = subs(d, [q0; v0], [q0_val; v0_val]);

% scale = double(subs(d, x, 0));
% n = n / scale;
% d = d / scale;
[n, d] = scaleQuotient(n, d);

% a = min(x0, 0);
% b = max(x0, 0);
% [prog, s] = prog.newSOSPoly(monomials(x, 0 : deg(n, x)));
% [prog, t] = prog.newSOSPoly(monomials(x, 0 : deg(n, x) - 2));
% prog = prog.withPolyEqs(s + (x - a) * (b - x) * t - n);

g_x = -x * (x - x0);

[prog, u_sos_1] = spotless_add_sprocedure(prog, d, g_x, x, even_degree(d, x) - even_degree(g_x, x));
[prog, u_sos_1] = spotless_add_sprocedure(prog, u_sos_1, n, x, even_degree(d, x) - even_degree(n, x));
prog = prog.withSOS(u_sos_1);

[prog, u_sos_2] = spotless_add_sprocedure(prog, -d, g_x, x, even_degree(d, x) - even_degree(g_x, x));
[prog, u_sos_2] = spotless_add_sprocedure(prog, u_sos_2, -n, x, even_degree(d, x) - even_degree(n, x));
prog = prog.withSOS(u_sos_2);

if ~isinf(u_max)
  u_max_sos = u_max * d - n;
  [prog, u_max_sos] = spotless_add_sprocedure(prog, u_max_sos, g_x, x, even_degree(u_max_sos, x) - even_degree(g_x, x));
  prog = prog.withSOS(u_max_sos);
end

cost = 0;

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(cost, solver, spot_options);

% ret = sol.isPrimalFeasible;
ret = sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
end

function [B, x, nu, du] = captureHeightControlSOS(xdot_sym, u_sym, q_sym, v_sym, g_val, zf_val)
checkDependency('spotless');
checkDependency('mosek');

prog = spotsosprog;
[prog, x] = prog.newIndeterminate('x', length(q_sym) + length(v_sym));
v = x(length(q_sym) + 1 : end);
[prog, vd] = prog.newIndeterminate('vd', length(v_sym));
xd = [v; vd];

[nu_sym, du_sym] = numden(u_sym);
du_sym = du_sym / q_sym(1)^4; 
nu = sym2msspoly([q_sym; v_sym], x, nu_sym);
du = sym2msspoly([q_sym; v_sym], x, du_sym);
[nu, du] = scaleQuotient(nu, du);

nvd = zeros(length(v_sym), 1) * msspoly(0);
dvd = zeros(length(v_sym), 1) * msspoly(0);
for i = 1 : length(v_sym)
  [nf_sym_i, df_sym_i] = numden(xdot_sym(i + length(q_sym)));
  nvd(i) = sym2msspoly([q_sym; v_sym], x, nf_sym_i);
  dvd(i) = sym2msspoly([q_sym; v_sym], x, df_sym_i);
  [nvd(i), dvd(i)] = scaleQuotient(nvd(i), dvd(i));
end
% dvd(1) = dvd(1) * x(1); nvd(1) = nvd(1) * x(1);

B_degree = 12;
epsilon = 1e-5;
[prog, B] = prog.newFreePoly(monomials(x, 0 : B_degree));

% B >= epsilon whenever nu >= 0 and du <= 0
B_unsafe_sos_1 = B - epsilon;
[prog, B_unsafe_sos_1] = spotless_add_sprocedure(prog, B_unsafe_sos_1, nu, x, B_degree - even_degree(nu, x));
[prog, B_unsafe_sos_1] = spotless_add_sprocedure(prog, B_unsafe_sos_1, -du, x, B_degree - even_degree(du, x));
prog = prog.withSOS(B_unsafe_sos_1);

% B >= epsilon whenever nu <= 0 and du >= 0
B_unsafe_sos_2 = B - epsilon;
[prog, B_unsafe_sos_2] = spotless_add_sprocedure(prog, B_unsafe_sos_2, -nu, x, B_degree - even_degree(nu, x));
[prog, B_unsafe_sos_2] = spotless_add_sprocedure(prog, B_unsafe_sos_2, du, x, B_degree - even_degree(du, x));
prog = prog.withSOS(B_unsafe_sos_2);

% TODO:
omega = sqrt(g_val / zf_val);
x0 = 1 * [-1; 1; omega; 0];
prog = prog.withPolyEqs(subs(B, x, x0) + 1);

% g_X0 = 1 - (x - x0)' * (x - x0) / 0.01^2;
% [prog, X0_sos] = spotless_add_sprocedure(prog, -B, g_X0, x, B_degree - even_degree(g_X0, x));
% prog = prog.withSOS(X0_sos);

Bdot_sos = -diff(B, x) * xd;
for i = 1 : length(vd)
  dynamics_eq = dvd(i) * vd(i) - nvd(i);
  [prog, Bdot_sos] = spotless_add_eq_sprocedure(prog, Bdot_sos, dynamics_eq, [x; vd], B_degree - even_degree(dynamics_eq, [x; vd]));
end
prog = prog.withSOS(Bdot_sos);

cost = 0;
% [prog, W] = prog.newSOSPoly(monomials(x, 0 : B_degree));
% prog = prog.withSOS(W - B - 1);
% 
% R_diag = 2 * ones(length(x), 1)';
% cost = spotlessIntegral(prog, W, x, R_diag, [], []);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(cost, solver, spot_options);
B = sol.eval(B);
end

function [xs, xds, feasible] = evaluateTrajectory(f, u_traj, u, xdot, q, v, q0, v0, q0_val, v0_val, g_val)
u_traj_val = subs(u_traj, [q0; v0], [q0_val; v0_val]);
feasible = captureHeightTrajectorySOS(u_traj_val, q(1), q0_val(1), inf);
disp(['feasible: ' num2str(feasible)])

f_val = subs(f, [q0; v0], [q0_val; v0_val]);
xs = linspace(q0_val(1), 0, 100);
zs = polyval(sym2poly(f_val), xs);
xds = horizontalVelocities(f_val, g_val, v0_val(1), xs);
dzdxs = polyval(sym2poly(diff(f_val, q(1))), xs);
zds = dzdxs .* xds;
f_legs = legForce(u, q, v, xs, zs, xds, zds);

half_index = floor(length(xs) / 2);
f_half_val = subs(f, [q0; v0], [xs(half_index); zs(half_index); xds(half_index); zds(half_index)]);
valuecheck(sym2poly(f_val), sym2poly(f_half_val));

xdot_fun = matlabFunction(xdot, 'Vars', [q; v]);
omega0 = sqrt(g_val / q0_val(2));
T = 3 / omega0;
x_init = [q0_val; v0_val];
[~, xtraj] = ode45(@(t, x) xdot_fun(x(1), x(2), x(3), x(4)), [0, T], x_init);

figure(1);
clf;

subplot(3, 1, 1);
hold on;
plot(xs, zs, 'b');
plot(xtraj(:, 1), xtraj(:, 2), 'r');
xlabel('q_x'); ylabel('q_z');
legend({'trajectory', 'simulation'}, 'Location', 'Best');
% axis(1.2 * [min(x0, 0), max(x0, 0), 0, max([z'; xtraj(:, 2)])]);
hold off;

subplot(3, 1, 2);
plot(xs, f_legs);
xlabel('q_x'); ylabel('f_{leg}');

subplot(3, 1, 3);
plot(xs, xds);
xlabel('q_x'); ylabel('v_x');
end

function [num, den] = scaleQuotient(num, den)
[~, ~, coefs] = decomp(den);
vd_scale = max(max(abs(coefs)));
num = num / vd_scale;
den = den / vd_scale;
end

function spot_poly = sym2msspoly(sym_vars, spot_vars, sym_poly)
[coefficients, monomials] = coeffs(sym_poly, sym_vars);
coefficients = double(coefficients);
spot_poly = msspoly(0);
for i = 1 : length(coefficients)
  term = msspoly(coefficients(i));
  for j = 1 : length(sym_vars)
    degree = double(feval(symengine, 'degree', monomials(i), sym_vars(j)));
    term = term * spot_vars(j)^degree;
  end
  spot_poly = spot_poly + term;
end
end
