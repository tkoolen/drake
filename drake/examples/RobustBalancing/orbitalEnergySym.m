function orbitalEnergySym()
n = 3;
syms g real;
syms x real;
syms xd real;
a = sym('a', [n + 1, 1]);
assumeAlso(a, 'real');

mons = sym.zeros(n + 1, 1);
for i = 0 : n
  mons(i + 1) = x^i;
end

f = a' * mons;
fp = diff(f, x); ap = coeffs(fp, x)';
fpp = diff(fp, x); app = coeffs(fpp, x)';
h = f - fp * x; ah = coeffs(h, x)';
fint = int(f * x, x); aint = coeffs(fint, x)';

s = diag(eye(n + 1));
sp = diag(double(jacobian([0; ap], a)));
spp = diag(double(jacobian([0; 0; app], a)));
sint = diag(double(jacobian(aint, a)));
sh = diag(double(jacobian([ah(1); 0; ah(2:end)], a)));

orbital_energy = simplify(1/2 * xd^2 * h^2 + g * x^2 * f - 3 * g * fint);
S = sh * sh';
r = g * (s - 3 * sint);

P = simplify(xd^2 * (sh .* mons) * (sh .* mons)'); %S .* (mons * mons'); % diff(diff(orbital_energy, a)', a);
q = simplify(x^2 * r .* mons); %diff(orbital_energy - 1/2 * a' * P * a, a)';
assert(logical(simplify(1/2 * a' * P * a + q' * a - orbital_energy) == 0));

% for E_orbit = 0:
force_numerator_divided_by_g = h^2 + 2 * fpp * (-x^2 * f + 3 * fint); % quadratic in coeffs of a; want this to be >= 0 on [x0, 0]
Pf = simplify(jacobian(jacobian(force_numerator_divided_by_g, a)', a));
qf = simplify(jacobian(force_numerator_divided_by_g - 1/2 * a' * Pf * a, a)');
valuecheck(double(qf), zeros(size(qf)));
Sf = s * s' + sp * sp' - s * sp' - sp * s' - s * spp' - spp * s' + 3 * spp * sint' + 3 * sint * spp';
assert(logical(simplify(mons' * (Sf .* (a * a')) * mons - force_numerator_divided_by_g) == 0));
% having M be positive semidefinite requires that a_n = 0 for n >= 3; if n < 3
% then it is trivially positive semidefinite

syms x0 xd0 real;
syms z0 zd0 zf real;
dzdx0 = zd0 / xd0;

C = [...
  jacobian(subs(f, x, x0), a); 
  jacobian(subs(fp, x, x0), a);
  jacobian(subs(f, x, 0), a)];
a_particular = C \ [z0; dzdx0; zf];
N = null(C);
syms w real;
orbital_energy_w = simplify(subs(orbital_energy, [x; xd; a], [x0; xd0; a_particular + N * w]));
[R, m] = coeffs(orbital_energy_w, w);
w = solve(orbital_energy_w == 0, w); % just a linear equation
a_sol = simplify(a_particular + N * w);
f_sol = simplify(subs(f, a, a_sol));
u = (g + fpp * xd^2) / (f - fp * x);
u_sol = simplify(subs(u, [x; xd; a], [x0; xd0; a_sol]));

g_val = 9.8;
x0_val = -1;
z0_val = 0.95;
omega0 = sqrt(g_val / z0_val);
xd0_val = -omega0 * x0_val * 0.95; %95; %1.05;
zd0_val = 0;
zf_val = 1;

a_val = double(subs(a_sol, [g; x0; z0; xd0; zd0; zf], [g_val; x0_val; z0_val; xd0_val; zd0_val; zf_val]));


end
