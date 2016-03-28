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

% find trajectory
[fw, u_trajw, uw, w, uhsq_trajw] = captureHeightTrajectory(g, q, v, q0, v0, zf, n);

nu = numden(uw);
syms a real;
nu_a = numden(simplify(subs(uw, x, a * xd)));
% nu_a = simplify(nu_a / x^4);
% nu_a_discriminant = evalin(symengine, ['polylib::discrim(' char(nu_a) ', a)']);
% nu_traj = simplify(subs(numden(u_trajw), xd0, a * x0));

% assumeAlso(g > 0);
% assumeAlso(zf > 0);
% a_roots = solve(nu_a == 0, a, 'MaxDegree', 4);

if n == 3
  syms c real;
  valuecheck(simplify(subs(fw, [x0; xd0; x], c * [x0; xd0; x]) - fw), sym(0));
  valuecheck(simplify(subs(uw, [x; xd], c * [x; xd]) - uw), sym(0));
  valuecheck(simplify(u_trajw - subs(u_trajw, [x0; xd0; x], c * [x0; xd0; x])), sym(0));
end

% plug in g and zf
g_val = 9.8;
zf_val = 1;
fw = subs(fw, [g; zf], [g_val; zf_val]);
uw = subs(uw, [g; zf], [g_val; zf_val]);
u_trajw = subs(u_trajw, [g; zf], [g_val; zf_val]);
uhsq_trajw = subs(uhsq_trajw, [g; zf], [g_val; zf_val]);
nu_a = subs(nu_a, [g; zf], [g_val; zf_val]);
% nu_a_discriminant = subs(nu_a_discriminant, [g; zf], [g_val; zf_val]);

% plug in w
w_val = randn(size(w));
f = subs(fw, w, w_val);
u = subs(uw, w, w_val);
u_traj = subs(u_trajw, w, w_val);
uhsq_traj = subs(uhsq_trajw, w, w_val);

% closed loop dynamics
vdot = [0; -g_val] + u * q;
xdot = simplify([v; vdot]);
% xdot = - 1 * [q;v];

% barrier function
% [B, state, nu, du] = captureHeightControlSOS(xdot, u, q, v, g_val, zf_val);

% gui setup
h_fig = orbitalEnergyGUI();
data = guidata(h_fig);
ax = data.axes1;
handles = guihandles(h_fig);
axes_callback = @(src, eventdata) plotTrajectory(src, eventdata, handles, f, u_traj, uhsq_traj, u, xdot, q, v, q0, v0, g_val);
set(ax,'ButtonDownFcn', axes_callback);
font_size = 16;
set(ax, 'FontSize', font_size);
set(findall(h_fig, 'type', 'Text'), 'FontSize', font_size);
slider_callback = @() plotStateSpaceSliceAndSetButtonDowns(handles, ax, u, u_traj, q, v, q0, v0, g_val, zf_val, w_val);
handles.z0slider.Callback = @(hObject, eventdata) slider_callback();
handles.zd0slider.Callback = @(hObject, eventdata) slider_callback();
end

function plotStateSpaceSliceAndSetButtonDowns(handles, ax, u, u_traj, q, v, q0, v0, g_val, zf_val, w_val)
plotStateSpaceSlice(ax, u, u_traj, q, v, q0, v0, g_val, zf_val, handles.z0slider.Value, handles.zd0slider.Value, w_val);
h = get(ax, 'children');
set(h, 'HitTest', 'off');
end

function plotTrajectory(src, eventdata, handles, f, u_traj, uhsq_traj, u, xdot, q, v, q0, v0, g_val)
x0 = eventdata.IntersectionPoint(1);
xd0 = eventdata.IntersectionPoint(2);
z0 = handles.z0slider.Value;
zd0 = handles.zd0slider.Value;
q0_val = [x0; z0];
v0_val = [xd0; zd0];

[xs, xds, feasible] = evaluateTrajectory(f, u_traj, u, xdot, q, v, q0, v0, q0_val, v0_val, g_val);

while ~isa(src, 'matlab.graphics.axis.Axes')
  src = src.Parent;
end
axes(src);
hold on;
if feasible
  color = 'g';
else
  color = 'r';
end
plot(xs(1), xds(1), [color 'x']);

% [verified, g_X0, x_xd_spot] = captureHeightTrajectoryRegionSOS(u_traj, q(1), [q0; v0], q0_val(1), q0_val(2), v0_val(1), v0_val(2), g_val);
% if verified
%   color = 'g';
% else
%   color = 'r';
% end
% V = axis(src);
% x_range = V(1 : 2);
% xd_range = V(3 : 4);
% contourSpotless(g_X0, x_xd_spot(1), x_xd_spot(2), x_range, xd_range, [], [], 0, {color});
% hold off;

end

function plotStateSpaceSlice(ax, u, u_traj, q, v, q0, v0, g_val, zf_val, z0_val, zd0_val, w_val)


cla(ax);
axes(ax);
hold on;

x = q(1);
z = q(2);
xd = v(1);
zd = v(2);

omega = sqrt(g_val / zf_val);
x_range = [-1, 1];
xd_range = x_range * omega;
axis([x_range, xd_range]);
axis manual
xs = linspace(x_range(1), x_range(2), 500);
xds = linspace(xd_range(1), xd_range(2), 500);

[Xs, Xds] = meshgrid(xs, xds);
u_fun = matlabFunction(subs(u, [z; zd],[z0_val; zd0_val]), 'Vars', [x; xd]);
Us = -u_fun(Xs, Xds);
[~, h_u] = contourf(Xs, Xds, Us, [0 0]);
valids = validTrajectory(g_val, zf_val, Xs, repmat(z0_val, size(Xs)), Xds, repmat(zd0_val, size(Xs)));
if any(any(valids))
  [~, h_valid] = contourf(Xs, Xds, valids, [0.5 inf]);
end

% contourSpotless(subs(B, [state(2); state(4)], [z0_val; zd0_val]), state(1), state(3), x_range, xd_range, [], [], 0, {'b'});
% contourSpotless(subs(-nu * du, [state(2); state(4)], [z0_val; zd0_val]), state(1), state(3), x_range, xd_range, [], [], 0, {'b'});
% contourSpotless(subs(du, [state(2); state(4)], [z0_val; zd0_val]), state(1), state(3), x_range, xd_range, [], [], 0, {'r'});
icp_line = -omega * xs;
h_icp = plot(xs, icp_line, 'b--');

% if all(w_val == 0)
%   if zd0_val == 0
%     a_crit = (-1/2).*5.^(-1/2).*(g_val.*z0_val.^(-2).*(7.*z0_val+3.*zf_val+3.^(1/2).*(3.*z0_val.^2+14.*z0_val.*zf_val+3.*zf_val.^2).^(1/2))).^(1/2);
%   else
%     [nu, du] = numden(u);
%     syms a real;
%     nu_a = simplify(subs(nu, xd, a * x));
%     nu_a = nu_a / (a * x^5);
%     nu_a = subs(nu_a, [z, zd], [z0_val, zd0_val]);
%     a_sols = roots(sym2poly(nu_a));
%     a_crit = min(a_sols);
%   end
%   plot(xs, a_crit * xs, 'r');
% %   plot(xs, max(a_sols) * xs, 'g');
% end

% [verified, g_X0, x0] = captureHeightTrajectoryRegionSOS(u_traj, q(1), [q0; v0], z0_val, zd0_val, g_val);
% if verified
%   contourSpotless(g_X0, x0(1), x0(2), x_range, xd_range, [], [], 0, {'g'});
% end

xlabel('$$x$$', 'Interpreter', 'latex');
ylabel('$$\dot{x}$$', 'Interpreter', 'latex');
legend([h_u, h_icp], {'$$u = 0$$', '$$x + \sqrt{\frac{z_f}{g}} \dot{x} = 0$$'}, 'Location', 'NorthEast', 'Interpreter', 'latex')
end

function [height_traj, u_traj, controller, w, uhsq] = captureHeightTrajectory(g, q, v, q0, v0, zf, n)
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

a_sol = A \ b + null(A) * w; % A \ b and null(A) seem to be computed using rref; A \ b is not a minimum norm solution!

h = f - fp * x;
fint = int(f * x, x);
% orbital_energy = simplify(1/2 * xd^2 * h^2 + g * x^2 * f - 3 * g * fint);

height_traj = simplify(subs(f, a, a_sol));
% u = (g + fpp * xd^2) / (z - fp * x);
u = (g + fpp * xd^2) / (z - fp * x);
u = simplify(subs(u, a, a_sol));
controller = simplify(subs(u, [x0; z0; xd0; zd0], [x; z; xd; zd]));

xd_squared = simplify(2 * (3 * g * fint - g * x^2 * f) / h^2);
u_traj = simplify((g + fpp * xd_squared) / (f - fp * x));
u_traj = simplify(subs(u_traj, a, a_sol));

uhsq = simplify(subs(simplify(g * h + fpp * g * 2 * int(h * x, x) / h), a, a_sol));

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

[n, d] = scaleQuotient(n, d);

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

ret = sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
end

function [verified, g_X0, x0, x0_sym] = captureHeightTrajectoryRegionSOS(utraj_sym, q_x_sym, x0_sym, x0_val, z0_val, xd0_val, zd0_val, g_val)
checkDependency('spotless');
checkDependency('mosek');

% utraj_sym = simplify(subs(utraj_sym, [x0_sym(2), x0_sym(4)], [z0_val, zd0_val]));
% x0_sym = [x0_sym(1); x0_sym(3)];
xc = [xd0_val / x0_val; z0_val; zd0_val];

% utraj_sym = simplify(partfrac(subs(utraj_sym, x0_sym(1), a * x0_sym(3))));
% a_val = xd0_val / x0_val;
% x0_sym = [x0_sym(1); a];
% xc = [x0_val; a_val];
% xc = [x0_val; z0_val; xd0_val; zd0_val];
syms a b real;
utraj_sym = simplify(subs(utraj_sym, [x0_sym(1); x0_sym(3); q_x_sym], [1; a; b]));
x0_sym(3) = a;
x0_sym(1) = [];
q_x_sym = b;

prog = spotsosprog;
[prog, x] = prog.newIndeterminate('x');
[prog, x0] = prog.newIndeterminate('xi', length(x0_sym));
% [prog, x0] = prog.newIndeterminate('xi', length(x0_sym));


[n_sym, d_sym] = numden(utraj_sym);
% n_sym = simplify(n_sym / x0_sym(1)^4); % NOTE: IMPORTANT when using UTRAJ

n = sym2msspoly([q_x_sym; x0_sym], [x; x0], n_sym);
d = sym2msspoly([q_x_sym; x0_sym], [x; x0], d_sym);

% n = subs(n, x0(3), zd0_val);
% d = subs(d, x0(3), zd0_val);
% x0 = [x0(1); x0(2)];

g_X = x * (1 - x);

% [~, ~, coefs] = decomp(n);
% scale = max(max(abs(coefs)));
% n = n / scale;

% [n, d] = scaleQuotient(n, d);

% sos_degree = 16;
% % g_X0_degree = 2;
% % [prog, g_X0] = prog.newFreePoly(monomials(delta, 0 : g_X0_degree));
% [prog, R] = prog.newFree(1);
% g_X0 = R - delta' * delta;
% % g_X0 = 0.2^2 - delta' * delta;
% 
% 
% % n negative, d positive, and g in [x0, 0] implies g_X0 < epsilon
% epsilon = 0; %1e-5;
% g_X0_sos_1 = -epsilon - g_X0;
% [prog, g_X0_sos_1] = spotless_add_sprocedure(prog, g_X0_sos_1, -n, [x; delta], sos_degree - even_degree(n, [x; delta]));
% [prog, g_X0_sos_1] = spotless_add_sprocedure(prog, g_X0_sos_1, g_X, [x; delta], sos_degree - even_degree(g_X, [x; delta]));
% [prog, g_X0_sos_1] = spotless_add_sprocedure(prog, g_X0_sos_1, 1 - [x; delta]' * [x; delta], [x; delta], sos_degree - 2);
% [prog, g_X0_sos_1] = spotless_add_sprocedure(prog, g_X0_sos_1, d, [x; delta], sos_degree - even_degree(d, [x; delta]));
% prog = prog.withSOS(g_X0_sos_1);
% 
% % g_X0_sos_2 = g_X0;
% % [prog, g_X0_sos_2] = spotless_add_sprocedure(prog, g_X0_sos_2, g_X, [x; x0], sos_degree - even_degree(n, [x; x0]));
% % [prog, g_X0_sos_2] = spotless_add_sprocedure(prog, g_X0_sos_2, -n, [x; x0], sos_degree - even_degree(n, [x; x0]));
% % [prog, g_X0_sos_2] = spotless_add_sprocedure(prog, g_X0_sos_2, d, [x; x0], sos_degree - even_degree(d, [x; x0]));
% % prog = prog.withSOS(g_X0_sos_2);
% 
% cost = -R;

% cost = subs(g_X0, x0, xc);
% R_diag = ones(size(delta))';
% [prog, w] = prog.newSOSPoly(monomials(delta, 0 : g_X0_degree));
% prog = prog.withSOS(w - g_X0 - 1);
% cost = spotlessIntegral(prog, w, delta, R_diag, [], []);

% works:
cost = 0;
radius = 0.1;
delta = x0 - xc;
g_X0 = radius^2 - delta' * delta;

% scale = abs(double(subs(n, [x0; x], [xc; xc(1)])));
% n = n / scale;
n_sos = n;
[prog, n_sos] = spotless_add_sprocedure(prog, n_sos, g_X0, [x; x0], even_degree(n_sos, [x; x0]) - even_degree(g_X0, x0));
[prog, n_sos] = spotless_add_sprocedure(prog, n_sos, g_X, [x; x0], even_degree(n_sos, [x; x0]) - even_degree(g_X, x0));
prog = prog.withSOS(n_sos);

% d_sos = d;
% [prog, d_sos] = spotless_add_sprocedure(prog, d_sos, g_X0, [x; x0], even_degree(n_sos, [x; x0]) - even_degree(g_X0, x0));
% [prog, d_sos] = spotless_add_sprocedure(prog, d_sos, g_X, [x; x0], even_degree(n_sos, [x; x0]) - even_degree(g_X, x0));
% prog = prog.withSOS(d_sos);

% % x and xd should have opposite sign:
% x_xd_sos = -g_X0;
% sos_degree = 10;
% [prog, x_xd_sos] = spotless_add_sprocedure(prog, x_xd_sos, prod(xc + delta), delta, sos_degree - 2);
% prog = prog.withSOS(x_xd_sos);

% solve
spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(cost, solver, spot_options);

verified = sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
if verified
  g_X0 = clean(sol.eval(g_X0));
end
% g_X0 = subs(g_X0, delta, x0 - xc);


end

function [B, x, nu, du] = captureHeightControlSOS(xdot_sym, u_sym, q_sym, v_sym, g_val, zf_val)
checkDependency('spotless');
checkDependency('mosek');


% reparameterize to have the origin correspond to z = zf:
syms dz real;
xdot_sym = subs(xdot_sym, q_sym(2), zf_val + dz);
u_sym = subs(u_sym, q_sym(2), zf_val + dz);
q_sym = [q_sym(1); dz];

xdot_sym = xdot_sym - 10 * [q_sym; v_sym];


prog = spotsosprog;
[prog, x] = prog.newIndeterminate('x', length(q_sym) + length(v_sym));
q = x(1 : length(q_sym));
v = x(length(q_sym) + 1 : end);

% state constraint
R_diag = 1 * ones(length(x), 1)';
A = diag(1./(R_diag.^2));
g_X = 1 - x'*A*x;

[nu_sym, du_sym] = numden(u_sym);
% du_sym = du_sym / q_sym(1)^4; 
nu = sym2msspoly([q_sym; v_sym], x, nu_sym);
du = sym2msspoly([q_sym; v_sym], x, du_sym);
[nu, du] = scaleQuotient(nu, du);

% xf = [0.5; 2; 0.5; 0];
% nu = 1 - (x - xf)' * (x - xf) / 0.1^2;
% du = 1;

nvd = zeros(length(v_sym), 1) * msspoly(0);
dvd = zeros(length(v_sym), 1) * msspoly(0);
for i = 1 : length(v_sym)
  [nvd_sym_i, dvd_sym_i] = numden(xdot_sym(i + length(q_sym)));
  nvd(i) = sym2msspoly([q_sym; v_sym], x, nvd_sym_i);
  dvd(i) = sym2msspoly([q_sym; v_sym], x, dvd_sym_i);
  [nvd(i), dvd(i)] = scaleQuotient(nvd(i), dvd(i));
end

% B_degree = 4;
sos_degree = 12;
% epsilon = 1e-5;
% [prog, B] = prog.newFreePoly(monomials(x, 0 : B_degree));
% B = x' * x - 0.01;
omega = sqrt(g_val / zf_val);
x0 = 1 * [-0.5; 0; 0.5 * omega; 0];
B = (x - x0)' * (x - x0) / 0.1^2 - 1;


% % B >= epsilon whenever nu >= 0 and du <= 0
% B_unsafe_sos_1 = B - epsilon;
% [prog, B_unsafe_sos_1] = spotless_add_sprocedure(prog, B_unsafe_sos_1, nu, x, sos_degree - even_degree(nu, x));
% [prog, B_unsafe_sos_1] = spotless_add_sprocedure(prog, B_unsafe_sos_1, -du, x, sos_degree - even_degree(du, x));
% prog = prog.withSOS(B_unsafe_sos_1);
% 
% % B >= epsilon whenever nu <= 0 and du >= 0
% B_unsafe_sos_2 = B - epsilon;
% [prog, B_unsafe_sos_2] = spotless_add_sprocedure(prog, B_unsafe_sos_2, -nu, x, sos_degree - even_degree(nu, x));
% [prog, B_unsafe_sos_2] = spotless_add_sprocedure(prog, B_unsafe_sos_2, du, x, sos_degree - even_degree(du, x));
% prog = prog.withSOS(B_unsafe_sos_2);

% % B >= epsilon outside of ball
% B_X = B - epsilon;
% [prog, B_X] = spotless_add_sprocedure(prog, B_X, -g_X, x, sos_degree - even_degree(g_X, x));
% prog = prog.withSOS(B_X);
% 
% % B >= epsilon when denominator is positive
% B_dvd = B - epsilon;
% [prog, B_dvd] = spotless_add_sprocedure(prog, B_dvd, dvd, x, sos_degree - even_degree(dvd, x));
% prog = prog.withSOS(B_dvd);

% Initial condition:
% omega = sqrt(g_val / zf_val);
% x0 = 1 * [-0.5; 0; 0.5 * omega; 0];
% prog = prog.withPolyEqs(subs(B, x, x0) + 1);

% g_X0 = 1 - (x - x0)' * (x - x0) / 0.01^2;
% [prog, X0_sos] = spotless_add_sprocedure(prog, -B, g_X0, x, sos_degree - even_degree(g_X0, x));
% prog = prog.withSOS(X0_sos);

% Bdot <=0
% vd as indeterminate version
% [prog, vd] = prog.newIndeterminate('vd', length(v_sym));
% xd = [v; vd];
% Bdot_sos = -diff(B, x) * xd;
% for i = 1 : length(vd)
%   dynamics_eq = dvd(i) * vd(i) - nvd(i);
%   [prog, Bdot_sos] = spotless_add_eq_sprocedure(prog, Bdot_sos, dynamics_eq, [x; vd], sos_degree - even_degree(dynamics_eq, [x; vd]));
% end
% prog = prog.withSOS(Bdot_sos);

% clear denominator version
dvd(1) = dvd(1) * x(1); nvd(1) = nvd(1) * x(1);
valuecheck(double(dvd(1) - dvd(2)), 0)
dvd = dvd(1);
Bdot_sos = diff(B, q) * v * dvd - diff(B, v) * nvd; % assuming dvd is negative


% [prog, Bdot_sos1] = spotless_add_sprocedure(prog, Bdot_sos, dvd, x, even_degree(Bdot_sos, x) - even_degree(dvd, x));
% [prog, Bdot_sos1] = spotless_add_sprocedure(prog, Bdot_sos1, -g_X, x, sos_degree - even_degree(g_X, x));
% prog = prog.withSOS(Bdot_sos1);
% 
% [prog, Bdot_sos2] = spotless_add_sprocedure(prog, -Bdot_sos, -dvd, x, even_degree(Bdot_sos, x) - even_degree(dvd, x));
% [prog, Bdot_sos2] = spotless_add_sprocedure(prog, Bdot_sos2, -g_X, x, sos_degree - even_degree(g_X, x));
% prog = prog.withSOS(Bdot_sos2);

[prog, Bdot_sos1] = spotless_add_sprocedure(prog, Bdot_sos, -B, x, sos_degree - even_degree(g_X, x));
[prog, Bdot_sos1] = spotless_add_sprocedure(prog, Bdot_sos1, -x(1) * x(3), x, sos_degree - 2);
prog = prog.withSOS(Bdot_sos1);

% [prog, Bdot_sos2] = spotless_add_sprocedure(prog, -Bdot_sos, -dvd, x, even_degree(Bdot_sos, x) - even_degree(dvd, x));
% [prog, Bdot_sos2] = spotless_add_sprocedure(prog, Bdot_sos2, -g_X, x, sos_degree - even_degree(g_X, x));
% prog = prog.withSOS(Bdot_sos2);

cost = 0;
% [prog, W] = prog.newSOSPoly(monomials(x, 0 : B_degree));
% prog = prog.withSOS(W - B - 1);
% cost = spotlessIntegral(prog, W, x, R_diag, [], []);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(cost, solver, spot_options);
B = sol.eval(B);
end

function [xs, xds, feasible] = evaluateTrajectory(f, u_traj, u, xdot, q, v, q0, v0, q0_val, v0_val, g_val)
% u_traj_val = subs(u_traj, [q0; v0], [q0_val; v0_val]);
% feasible = q0_val(1) * v0_val(1) < 0 && captureHeightTrajectorySOS(u_traj_val, q(1), q0_val(1), inf);
% disp(['feasible: ' num2str(feasible)])

f_val = subs(f, [q0; v0], [q0_val; v0_val]);
xs = linspace(q0_val(1), 0, 100);
zs = polyval(sym2poly(f_val), xs);
xds = horizontalVelocities(f_val, g_val, v0_val(1), xs);
dzdxs = polyval(sym2poly(diff(f_val, q(1))), xs);
zds = dzdxs .* xds;
f_legs = legForce(u, q, v, xs, zs, xds, zds);
nfun = matlabFunction(numden(u), 'Vars', [q; v]);
ns = nfun(xs, zs, xds, zds);
feasible = q0_val(1) * v0_val(1) < 0 && all(f_legs(1 : end - 1) > 0);

% half_index = floor(length(xs) / 2);
% f_half_val = subs(f, [q0; v0], [xs(half_index); zs(half_index); xds(half_index); zds(half_index)]);
% valuecheck(sym2poly(f_val), sym2poly(f_half_val), 1e-3);

xdot_fun = matlabFunction(xdot, 'Vars', [q; v]);
zf = double(subs(f_val, q(1), 0));
omega0 = sqrt(g_val / zf);
T = 3 / omega0;
x_init = [q0_val; v0_val];
[ttraj, xtraj] = ode45(@(t, x) xdot_fun(x(1), x(2), x(3), x(4)), [0, T], x_init);

h_fig = figure(1);
clf;

subplot(3, 1, 1);
hold on;
plot(xs, zs, 'b');
plot(xtraj(:, 1), xtraj(:, 2), 'r--');
xlabel('$$x$$', 'Interpreter', 'latex');
ylabel('$$z$$', 'Interpreter', 'latex');
legend({'trajectory', 'simulation'}, 'Location', 'Best');
% axis(1.2 * [min(x0, 0), max(x0, 0), 0, max([z'; xtraj(:, 2)])]);
hold off;

subplot(3, 1, 2);
plot(xs, xds);
xlabel('$$x$$', 'Interpreter', 'latex');
ylabel('$$\dot{x}$$', 'Interpreter', 'latex');

subplot(3, 1, 3);
hold on
plot(xs, f_legs);
% plot(xs, ns, 'r');
hold off
xlabel('$$x$$', 'Interpreter', 'latex');
ylabel('$$\frac{f_{l}}{m}$$', 'Interpreter', 'latex');

% subplot(4, 1, 4);
% hold on
% plot(xs, ns, 'r');
% hold off
% xlabel('$$x$$', 'Interpreter', 'latex');
% ylabel('$$n(x)$$', 'Interpreter', 'latex');

font_size = 12;
set(findall(h_fig, 'type', 'Axes'), 'FontSize', font_size);

end

function [num, den] = scaleQuotient(num, den)
[~, ~, coefs] = decomp(den);
scale = max(max(abs(coefs)));
num = num / scale;
den = den / scale;
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

function ret = validTrajectory(g, zf, x0, z0, xd0, zd0)
a = xd0 ./ x0;
y = zd0 - a .* z0;

ret = ((a)<(0))&((7.*a.^(-1).*g+20.*y)>=(3.^(1/2).*(g.*(3.*a.^(-2).*g+ ...
  40.*zf)).^(1/2)));
end
