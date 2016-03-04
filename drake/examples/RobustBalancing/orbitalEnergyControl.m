function orbitalEnergyControl()
checkDependency('spotless');
checkDependency('snopt');

n = 5;
%g = msspoly('g');
g = 9.8;
x = msspoly('x');
xd = msspoly('xd');
a = msspoly('a', n + 1);

mons = monomials(x, 0 : n);
f = a' * mons;
fp = diff(f, x); [ap, pp] = pdecomp(fp, x); ap = ap';
fpp = diff(fp, x); [app, ppp] = pdecomp(fpp, x); app = app';
h = f - fp * x; [ah, ph] = pdecomp(h, x); ah = ah';
fint = integral(f * x, x); [aint, pint] = pdecomp(fint, x); aint = aint';
orbital_energy = 1/2 * xd^2 * h^2 + g * x^2 * f - 3 * g * fint;
P = diff(diff(orbital_energy, a)', a);
q = diff(orbital_energy - 1/2 * a' * P * a, a)';
valuecheck(double(orbitalEnergy(a, P, q) - orbital_energy), 0);

% for E_orbit = 0:
force_numerator_divided_by_g = h^2 + 2 * fpp * (-x^2 * f + 3 * fint); % quadratic in coeffs of a; want this to be >= 0 on [x0, 0]
Pf = diff(diff(force_numerator_divided_by_g, a)', a);
qf = diff(force_numerator_divided_by_g - 1/2 * a' * Pf * a, a)';
valuecheck(double(qf), zeros(size(qf)));

Mf = a * a' + [0; ap] * [0; ap]' - a * [0; ap]' - [0; ap] * a' - a * [0; 0; app]' - [0; 0; app] * a' + 3 * [0; 0; app] * aint' + 3 * aint * [0; 0; app]'; 
valuecheck(double(clean(mons' * Mf * mons - force_numerator_divided_by_g)), 0);
% having M be positive semidefinite requires that a_n = 0 for n >= 3; if n < 0
% then it is trivially positive semidefinite

x0 = -0.2;
z0 = 1;
omega0 = sqrt(g / z0);
xd0 = -omega0 * x0 * 1.1; %95; %1.05;
zd0 = 0;
dzdx0 = zd0 / xd0;
zf = 1;
P0 = double(subs(P, [x; xd], [x0; xd0]));
q0 = double(subs(q, [x; xd], [x0; xd0]));

prog = NonlinearProgram(n + 1);

orbital_energy_constraint = FunctionHandleConstraint(0, 0, n + 1, @(a) orbitalEnergy(a, P0, q0), 2);
prog = prog.addNonlinearConstraint(orbital_energy_constraint);

for x_knot = linspace(x0, 0, 2 * n) % just enforce at a couple of sample points
  Pf_knot = double(subs(Pf, x, x_knot));
  force_constraint = FunctionHandleConstraint(0, inf, n + 1, @(a) forceConstraint(a, Pf_knot), 2); % only require that initial force >= 0...
  prog = prog.addNonlinearConstraint(force_constraint);
end

initial_height_constraint = LinearConstraint(z0, z0, double(diff(subs(f, x, x0), a)));
prog = prog.addLinearConstraint(initial_height_constraint);

initial_slope_constraint = LinearConstraint(dzdx0, dzdx0, double(diff(subs(diff(f, x), x, x0), a)));
prog = prog.addLinearConstraint(initial_slope_constraint);

final_height_constraint = LinearConstraint(zf, zf, double(diff(subs(f, x, 0), a)));
prog = prog.addLinearConstraint(final_height_constraint);

prog = prog.addQuadraticCost(eye(n + 1), zeros(n + 1, 1));

a0 = randn(n + 1, 1);
a_sol = prog.solve(a0);

xs = linspace(x0, 0, 100);
zs = polyval(flipud(a_sol), xs);

zd0 = polyval(polyder(flipud(a_sol)), x0) * xd0;

T = 5 / omega0;
x_init = [x0; z0; xd0; zd0];
[ttraj, xtraj] = ode45(@(t, x) closedLoopDynamics(x, g, a_sol), [0, T], x_init);
utraj = zeros(length(ttraj), 1);
force_numerator_divided_by_g_traj = zeros(length(ttraj), 1);
for i = 1 : length(ttraj)
  utraj(i) = control(xtraj(i, :)', g, a_sol);
  force_numerator_divided_by_g_traj(i) = msubs(force_numerator_divided_by_g, [a; x], [a_sol; xtraj(i, 1)]);
end

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
plot(ttraj, sqrt(sum(xtraj(:, 3 : 4).^2, 2)));
xlabel('t'); ylabel('|v|')

subplot(3, 1, 3);
hold on;
plot(ttraj, utraj, 'r');
% plot(ttraj, force_numerator_divided_by_g_traj, 'b');
xlabel('t'); ylabel('u')

end

function [E, dE] = orbitalEnergy(a, P, q)
E = 1/2 * a' * P * a + q' * a;
if nargout > 1
  dE = a' * P + q';
end
end

function [v, dv] = forceConstraint(a, Pf)
v = 1/2 * a' * Pf * a;
if nargout > 1
  dv = a' * Pf;
end
end

function u = computeForce(traj_coeffs, g, q, xd_squared)
x = q(1);
z = q(2);
traj_coeffs = flipud(traj_coeffs);
fp = polyder(traj_coeffs);
fpp = polyder(fp);
u = norm(q) * (g + polyval(fpp, x) * xd_squared) / (z - polyval(fp, x) * x);
end

function u = control(state, g, traj_coeffs)
q = state(1:2);
v = state(3:4);
x = q(1);
z = q(2);
xd = v(1);
xd_squared = xd^2;
u = computeForce(traj_coeffs, g, q, xd_squared);

k = 10;
u = u + k * (polyval(flipud(traj_coeffs), x) - z);
end

function xd = dynamics(x, u, g)
q = x(1:2);
v = x(3:4);
vd = u * q / norm(q) - [0; g];
xd =  [v; vd];
end

function xd = closedLoopDynamics(x, g, f)
u = control(x, g, f);
xd = dynamics(x, u, g);
end

