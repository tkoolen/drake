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

[f_sol, u_sol] = cubicCaptureHeightTrajectory(g, q, v, q0, v0, zf);

vdot_sol = [0; -g] + u_sol * q;
xdot_sol = simplify([v; vdot_sol]);

g_val = 9.8;
x0_val = -1;
z0_val = 0.95;
omega0 = sqrt(g_val / z0_val);
xd0_val = -omega0 * x0_val * 0.95; %1.05;
zd0_val = 0;
zf_val = 1;

f_val = subs(f_sol, [g; x0; z0; xd0; zd0; zf], [g_val; x0_val; z0_val; xd0_val; zd0_val; zf_val]);
u_val = subs(u_sol, [g; zf], [g_val; zf_val]);
xs = linspace(x0_val, 0, 100);
zs = polyval(sym2poly(f_val), xs);
xds = computeHorizontalVelocities(f_val, g_val, xd0_val, xs);
dzdxs = polyval(sym2poly(diff(f_val, x)), xs);
zds = dzdxs .* xds;
us = computeInputs(u_val, q, v, xs, zs, xds, zds);

half_index = floor(length(xs) / 2);
f_half_val = subs(f_sol, [g; x0; z0; xd0; zd0; zf], [g_val; xs(half_index); zs(half_index); xds(half_index); zds(half_index); zf_val]);
valuecheck(sym2poly(f_val), sym2poly(f_half_val));

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
legend({'trajectory', 'simulation'});
% axis(1.2 * [min(x0, 0), max(x0, 0), 0, max([z'; xtraj(:, 2)])]);
hold off;

subplot(3, 1, 2);
plot(xs, us);
xlabel('q_x'); ylabel('u');

subplot(3, 1, 3);
plot(xs, xds);
xlabel('q_x'); ylabel('v_x');

end

function [f_sol, u_sol] = cubicCaptureHeightTrajectory(g, q, v, q0, v0, zf)
x = q(1);
z = q(2);
xd = v(1);
zd = v(2);
x0 = q0(1);
z0 = q0(2);
xd0 = v0(1);
zd0 = v0(2);

n = 3;
a = sym('a', [n + 1, 1]);
assumeAlso(a, 'real');

mons = sym.zeros(n + 1, 1);
for i = 0 : n
  mons(i + 1) = x^i;
end

f = a' * mons;
fp = diff(f, x);
fpp = diff(fp, x);
h = f - fp * x;
fint = int(f * x, x);

orbital_energy = simplify(1/2 * xd^2 * h^2 + g * x^2 * f - 3 * g * fint);

dzdx0 = zd0 / xd0;

C = [...
  jacobian(subs(f, x, x0), a); 
  jacobian(subs(fp, x, x0), a);
  jacobian(subs(f, x, 0), a)];
a_particular = C \ [z0; dzdx0; zf];
N = null(C);
syms w real;
orbital_energy_w = simplify(subs(orbital_energy, [x; xd; a], [x0; xd0; a_particular + N * w]));
% [R, m] = coeffs(orbital_energy_w, w);
w = solve(orbital_energy_w == 0, w); % just a linear equation
a_sol = simplify(a_particular + N * w);
f_sol = simplify(subs(f, a, a_sol));

u = (g + fpp * xd^2) / (z - fp * x);
u_sol = simplify(subs(u, a, a_sol));
u_sol = simplify(subs(u_sol, [x0; z0; xd0; zd0], [x; z; xd; zd]));

% ap = coeffs(fp, x)';
% app = coeffs(fpp, x)';
% ah = coeffs(h, x)';
% aint = coeffs(fint, x)';
% 
% s = diag(eye(n + 1));
% sp = diag(double(jacobian([0; ap], a)));
% spp = diag(double(jacobian([0; 0; app], a)));
% sint = diag(double(jacobian(aint, a)));
% sh = diag(double(jacobian([ah(1); 0; ah(2:end)], a)));
% 
% S = sh * sh';
% r = g * (s - 3 * sint);
% 
% P = simplify(xd^2 * (sh .* mons) * (sh .* mons)'); %S .* (mons * mons'); % diff(diff(orbital_energy, a)', a);
% q = simplify(x^2 * r .* mons); %diff(orbital_energy - 1/2 * a' * P * a, a)';
% assert(logical(simplify(1/2 * a' * P * a + q' * a - orbital_energy) == 0));
% 
% % for E_orbit = 0:
% force_numerator_divided_by_g = h^2 + 2 * fpp * (-x^2 * f + 3 * fint); % quadratic in coeffs of a; want this to be >= 0 on [x0, 0]
% Pf = simplify(jacobian(jacobian(force_numerator_divided_by_g, a)', a));
% qf = simplify(jacobian(force_numerator_divided_by_g - 1/2 * a' * Pf * a, a)');
% valuecheck(double(qf), zeros(size(qf)));
% Sf = s * s' + sp * sp' - s * sp' - sp * s' - s * spp' - spp * s' + 3 * spp * sint' + 3 * sint * spp';
% assert(logical(simplify(mons' * (Sf .* (a * a')) * mons - force_numerator_divided_by_g) == 0));
% % having M be positive semidefinite requires that a_n = 0 for n >= 3; if n < 3
% % then it is trivially positive semidefinite
end

function xds = computeHorizontalVelocities(f, g, xd0, xs)
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

function us = computeInputs(u, q, v, xs, zs, xds, zds)
fun = matlabFunction(u, 'Vars', [q; v]);
us = fun(xs, zs, xds, zds);
end
