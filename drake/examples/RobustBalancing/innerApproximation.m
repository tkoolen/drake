function [Vsol,Wsol] = innerApproximation(model,u,R_diag,target,options)

if ~isfield(options,'beta')
  options.beta = 1;
end

%scaling of state vector
if ~isfield(options,'scale')
  scale = ones(model.num_states,1);
elseif length(options.scale) == 1
  scale = ones(model.num_states,1)*options.scale;
else
  scale = options.scale;
end
scale_inv = 1./scale;

degree = options.degree; % degree of V,W

%% Create SOS program
prog = spotsosprog;


%% Create indeterminates

[prog,x]=prog.newIndeterminate('x', model.num_states); % state

%% Scale r_diag
R_diag = scale'.*R_diag;

%% Create polynomials V(t,x) and W(x)
nbeta = length(options.beta);
V_vars = x;
[prog,V] = prog.newFreePoly(monomials(V_vars,0:degree),nbeta);

W_vars = x;
[prog,W] = prog.newFreePoly(monomials(W_vars,0:degree));

%% Dynamics
f = scale.*model.dynamics([], scale_inv.*x, u);

Vdot = diff(V,x)*f;

Vdot_degree = even_degree(Vdot,x);

%% Goal region

V0p = target(scale_inv.*x);

% State constraint
A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x;


%% SOS constraints

% (1) w(x) >= 1 + sum(V) for x not in goal region
[prog, wv_sos] = spotless_add_sprocedure(prog, W-sum(V) - 1, -V0p,V_vars,degree-2);
[prog, wv_sos] = spotless_add_sprocedure(prog, wv_sos, h_X,V_vars,degree-2);
prog = prog.withSOS(wv_sos);

% (2) w(x) >=0 for x not in goal region
[prog, w_sos] = spotless_add_sprocedure(prog, W, -V0p,V_vars,degree-2);
[prog, w_sos] = spotless_add_sprocedure(prog, w_sos, h_X,V_vars,degree-2);
prog = prog.withSOS(w_sos);

% (3) V(x) >=0 for x on boundary of state space
for i=1:nbeta,
  [prog, v_sos] = spotless_add_eq_sprocedure(prog, V(i), h_X,V_vars,degree-2);
  prog = prog.withSOS(v_sos);
end

% (4) Vdot(x) <= beta*v for x not in goal region
for i=1:nbeta,
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, options.beta(i)*V(i) - Vdot, -V0p,V_vars,Vdot_degree-2);
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, h_X,V_vars,Vdot_degree-2);
  prog = prog.withSOS(Vdot_sos);
end

%% Set up cost function -- integration over a sphere

cost = spotlessIntegral(prog,W,x,R_diag,[],[]);

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = false;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(cost,solver,spot_options);



%% Plotting
Vsol = subs(sol.eval(V),x,scale.*x);
Wsol = subs(sol.eval(W),x,scale.*x);
model.plotfun(0, sum(Vsol), Wsol, subs(h_X,x,scale.*x), R_diag, msspoly('t',1), x, []);
sum(Vsol)
end