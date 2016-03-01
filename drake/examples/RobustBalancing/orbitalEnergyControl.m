function orbitalEnergyControl()
checkDependency('spotless');
checkDependency('snopt')

n = 3;
g = 9.8;
x = msspoly('x');
xd = msspoly('xd');
a = msspoly('a', n + 1);
f = a' * monomials(x, 0 : n);
h = f - diff(f, x) * x;
orbital_energy = 1/2 * xd^2 * h^2 + g * x^2 * f - 3 * g * integral(f * x, x);
P = diff(diff(orbital_energy, a)', a);
q = diff(orbital_energy - 1/2 * a' * P * a, a)';
valuecheck(double(orbitalEnergy(a, P, q) - orbital_energy), 0)


x0 = -1;
z0 = 1;
omega0 = sqrt(g / z0);
xd0 = -omega0 * x0 + 0.1;
zd0 = 0;
dzdx0 = zd0 / xd0;
zf = 1;
P0 = double(subs(P, [x; xd], [x0; xd0]));
q0 = double(subs(q, [x; xd], [x0; xd0]));

prog = NonlinearProgram(n + 1);
orbital_energy_constraint = FunctionHandleConstraint(0, 0, n + 1, @(a) orbitalEnergy(a, P0, q0), 2);
prog = prog.addNonlinearConstraint(orbital_energy_constraint);

initial_height_constraint = LinearConstraint(z0, z0, double(diff(subs(f, x, x0), a)));
prog = prog.addLinearConstraint(initial_height_constraint);

initial_slope_constraint = LinearConstraint(dzdx0, dzdx0, double(diff(subs(diff(f, x), x, x0), a)));
prog = prog.addLinearConstraint(initial_slope_constraint);

final_height_constraint = LinearConstraint(zf, zf, double(diff(subs(f, x, 0), a)));
prog = prog.addLinearConstraint(final_height_constraint);

prog = prog.addQuadraticCost(eye(n + 1), zeros(n + 1, 1));

a0 = randn(n + 1, 1);
a_sol = prog.solve(a0);

x = linspace(x0, 0, 100);
z = polyval(flipud(a_sol), x);
plot(x, z, 'b');

T = 5 / omega0;
x_init = [x0; z0; xd0; zd0];
[~, xtraj] = ode45(@(t, x) closedLoopDynamics(x, g, a_sol), [0, T], x_init);
hold on;
plot(xtraj(:, 1), xtraj(:, 2), 'r');

end

function [E, dE] = orbitalEnergy(a, P, q)
E = 1/2 * a' * P * a + q' * a;
if nargout > 1
  dE = a' * P + q';
end
end

function u = control(state, g, f)
q = state(1:2);
v = state(3:4);
f = flipud(f);
x = q(1);
z = q(2);
xd = v(1);
fp = polyder(f);
fpp = polyder(fp);
u = norm(q) * (g + polyval(fpp, x) * xd^2) / (z - polyval(fp, x) * x);
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

