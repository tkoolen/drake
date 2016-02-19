function B_sol = barrierFunctionForClosedLoopSystem(model, R_diag, t, x, u, options)

xdot_closed_loop = model.dynamics(t, x, u);

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(t);
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
Bdot = diff(B, x) * xdot_closed_loop + diff(B, t);
prog = prog.withSOS(-Bdot);

% Solve
spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(subs(B, x, x0), solver, spot_options);
B_sol = sol.eval(B);

end
