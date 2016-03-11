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
n = 3;

g_val = 9.8;
zf_val = 1;

[fw, u_trajw, controllerw, w] = captureHeightTrajectory(g, q, v, q0, v0, zf, n);
u_trajw = subs(u_trajw, [g; zf], [g_val; zf_val]);
uw = subs(u_trajw, [q0; v0], [q; v]);

% w_val = captureHeightControlSOS(uw, q, v, w, inf);

x0_val = -1;
z0_val = 1; %0.95;
omega0 = sqrt(g_val / z0_val);
xd0_val = -omega0 * x0_val * 0.95; %1.05;
zd0_val = 0;
q0_val = [x0_val; z0_val];
v0_val = [xd0_val; zd0_val];
captureHeightTrajectorySOS(u_trajw, x, w, q0, v0, q0_val, v0_val, inf);
w_val = 0 * ones(size(w));
f_sol = subs(fw, w, w_val);
u_sol = subs(uw, w, w_val);

vdot_sol = [0; -g] + u_sol * q;
xdot_sol = simplify([v; vdot_sol]);
% [nxdot, dxdot] = numden(xdot_sol);




f_val = subs(f_sol, [g; x0; z0; xd0; zd0; zf], [g_val; x0_val; z0_val; xd0_val; zd0_val; zf_val]);
u_val = subs(u_sol, [g; zf], [g_val; zf_val]);
xs = linspace(x0_val, 0, 100);
zs = polyval(sym2poly(f_val), xs);
xds = horizontalVelocities(f_val, g_val, xd0_val, xs);
dzdxs = polyval(sym2poly(diff(f_val, x)), xs);
zds = dzdxs .* xds;
f_legs = legForce(u_val, q, v, xs, zs, xds, zds);

half_index = floor(length(xs) / 2);
f_half_val = subs(f_sol, [g; x0; z0; xd0; zd0; zf], [g_val; xs(half_index); zs(half_index); xds(half_index); zds(half_index); zf_val]);
% valuecheck(sym2poly(f_val), sym2poly(f_half_val));

xdot_fun = matlabFunction(subs(xdot_sol, [g; zf], [g_val; zf_val]), 'Vars', [q; v]);
T = 3 / omega0;
x_init = [x0_val; z0_val; xd0_val; zd0_val];
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
% u_traj_den = simplify(subs(h^2 * g + fpp * 2 * (3 * g * fint - g * x^2 * f), a, a_sol));
% u_traj_num = simplify(subs(h^3, a, a_sol));

u_traj = simplify((g + fpp * xd_squared) / (f - fp * x));
u_traj = simplify(subs(u_traj, a, a_sol));

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

function captureHeightTrajectorySOS(uwtraj_sym, x_sym, w_sym, q0_sym, v0_sym, q0_val, v0_val, u_max)
checkDependency('spotless');
checkDependency('mosek');

prog = spotsosprog;
[prog, x] = prog.newIndeterminate('x');
[prog, w] = prog.newFree(length(w_sym));

uwtraj_sym = subs(uwtraj_sym, [q0_sym; v0_sym], [q0_val; v0_val]);
[n_sym, d_sym] = numden(uwtraj_sym);
n = sym2msspoly([x_sym; w_sym], [x; w], n_sym);
d = sym2msspoly([x_sym; w_sym], [x; w], d_sym);

% q0 = msspoly('y', length(q0_sym));
% v0 = msspoly('z', length(q0_sym));
% n = sym2msspoly([x_sym; q0_sym; v0_sym; w_sym], [x; q0; v0; w], n_sym);
% d = sym2msspoly([x_sym; q0_sym; v0_sym; w_sym], [x; q0; v0; w], d_sym);
% n = subs(n, [q0; v0], [q0_val; v0_val]);
% d = subs(d, [q0; v0], [q0_val; v0_val]);

scale = double(subs(d, x, 0));
n = n / scale;
d = d / scale;

x0 = q0_val(1);
g_x = -x * (x - x0);

[prog, n_pos_sos] = spotless_add_sprocedure(prog, n, g_x, x, even_degree(n, x) - even_degree(g_x, x));
prog = prog.withSOS(n_pos_sos);

[prog, d_pos_sos] = spotless_add_sprocedure(prog, d, g_x, x, even_degree(d, x) - even_degree(g_x, x));
prog = prog.withSOS(d_pos_sos);

if ~isinf(u_max)
  u_max_sos = -(n - u_max * d);
  [prog, u_max_sos] = spotless_add_sprocedure(prog, u_max_sos, g_x, x, even_degree(n, x) - even_degree(g_x, x));
  prog = prog.withSOS(u_max_sos);
end

cost = 0;

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(cost, solver, spot_options);

% w = sol.eval(w);
end

function w = captureHeightControlSOS(uw_sym, q_sym, v_sym, w_sym, u_max)
checkDependency('spotless');
checkDependency('mosek');

prog = spotsosprog;
[prog, x] = prog.newIndeterminate('q', length(q_sym) + length(v_sym));
[prog, w] = prog.newFree(length(w_sym));

[nw_sym, d_sym] = numden(uw_sym);
nw = sym2msspoly([q_sym; v_sym; w_sym], [x; w], nw_sym);
d = sym2msspoly([q_sym; v_sym], x, d_sym);

% g_x = 1 - x' * x / 0.1^2; % TODO
% x0 = [-1; 1; 3.1305; 0];
% g_x = 1 - (x - x0)' * (x - x0) / 0.1^2;
omega = sqrt(9.8 / 1);
icp = x(1) + x(3) / omega;
g_x = [(0 * icp - icp^2); 1 - (x(2) - 1)^2 / 0.5^2; 1 - x(4)^2 / 0.5^2];

[prog, nw_pos_sos] = spotless_add_sprocedure(prog, nw, g_x, x, even_degree(nw, x) - even_degree(g_x, x));
[prog, nw_pos_sos] = spotless_add_sprocedure(prog, nw_pos_sos, d, x, even_degree(nw, x) - even_degree(d, x));
prog = prog.withSOS(nw_pos_sos);

[prog, nw_neg_sos] = spotless_add_sprocedure(prog, -nw, g_x, x, even_degree(nw, x) - even_degree(g_x, x));
[prog, nw_neg_sos] = spotless_add_sprocedure(prog, nw_neg_sos, -d, x, even_degree(nw, x) - even_degree(d, x));
prog = prog.withSOS(nw_neg_sos);

if ~isinf(u_max)
  u_max_sos = -(nw - u_max * d);
  [prog, u_max_sos] = spotless_add_sprocedure(prog, u_max_sos, g_x, x, even_degree(nw, x) - even_degree(g_x, x));
  prog = prog.withSOS(u_max_sos);
end

cost = 0;

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(cost, solver, spot_options);

w = sol.eval(w);

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
