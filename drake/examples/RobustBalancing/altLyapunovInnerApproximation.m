function [Vsol,Vlvl] = altLyapunovInnerApproximation(model,W,T,u,R_diag,options)

mode = 0;

%scaling of state vector
if ~isfield(options,'scale')
  scale = ones(model.num_states,1);
elseif length(options.scale) == 1
  scale = ones(model.num_states,1)*options.scale;
else
  scale = options.scale;
end

if ~isfield(options,'beta')
  beta = 0;
else
  beta = options.beta;
end

if ~isfield(options,'u_den')
  options.u_den = msspoly(1);
end

if ~isfield(options,'target')
  options.target = [];
end

if isfield(options,'Wlvl')
  Wlvl = options.Wlvl;
else
  Wlvl = 0;
end

target = options.target;
scale_inv = 1./scale;
u_den = options.u_den;
degree = options.degree; % degree of V,W

%% Create SOS program
prog = spotsosprog;


%% Create indeterminates
[prog,x]=prog.newIndeterminate('x', model.num_states); % state
[prog,t]=prog.newIndeterminate('t', 1); % time
u_mss = msspoly('u',model.num_inputs);

B = 1 - x'*x;

%% Scale r_diag


%% Create polynomials V(t,x) and W(x)
V_vars = [t;x];
[prog,V] = prog.newFreePoly(monomials(V_vars,1:degree));

% [prog,q] = prog.newFreePoly(monomials(V_vars,0:degree));

%% Dynamics
f = scale.*model.dynamics([], scale_inv.*x, u_mss);

T_init = T;
f = f*T;
T = 1;

g = diff(f,u_mss);

Vdot = diff(V,x)*subs(f,u_mss,u_mss*0)*u_den + diff(V,t)*u_den + diff(V,x)*g*u;

Vdot_degree = even_degree(Vdot,V_vars);
W_degree = even_degree(W,V_vars);
%% Goal region

% State constraint
R_diag = scale'.*R_diag;
A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x;


%% SOS constraints

% (1) w(x) >= Wlvl, x in X ==> Vdot <= 0
[prog, vdot_sos] = spotless_add_sprocedure(prog, -Vdot-beta*V, h_X,V_vars,Vdot_degree-2);
[prog, vdot_sos] = spotless_add_sprocedure(prog, vdot_sos, W-Wlvl,V_vars,Vdot_degree-W_degree);
[prog, vdot_sos] = spotless_add_sprocedure(prog, vdot_sos, t*(1-t),V_vars,Vdot_degree-2);
if ~isempty(target)
  [prog, vdot_sos] = spotless_add_sprocedure(prog, vdot_sos, -target(x),V_vars,Vdot_degree-2);
end
prog = prog.withSOS(vdot_sos);

if mode == 1

% (2) w(x) >= Wlvl, x in X ==> V <= 1
[prog, vmax_sos] = spotless_add_sprocedure(prog, 1-V, h_X,V_vars,degree-2);
[prog, vmax_sos] = spotless_add_sprocedure(prog, vmax_sos, W-Wlvl,V_vars,degree-W_degree);
[prog, vmax_sos] = spotless_add_sprocedure(prog, vmax_sos, t*(1-t),V_vars,degree-2);
prog = prog.withSOS(vmax_sos);

else
% (2) w(x) = Wlvl, x in X ==> V >= 1
[prog, vmax_sos] = spotless_add_sprocedure(prog, V - 1, h_X,V_vars,degree-2);
[prog, vmax_sos] = spotless_add_eq_sprocedure(prog, vmax_sos, W-Wlvl,V_vars,degree-W_degree);
[prog, vmax_sos] = spotless_add_sprocedure(prog, vmax_sos, t*(1-t),V_vars,degree-2);
prog = prog.withSOS(vmax_sos);

% (2b) w(x) >= Wlvl, x in dX ==> V >= 1
[prog, vmax_sos] = spotless_add_eq_sprocedure(prog, V - 1, h_X,V_vars,degree-2);
[prog, vmax_sos] = spotless_add_sprocedure(prog, vmax_sos, W-Wlvl,V_vars,degree-W_degree);
[prog, vmax_sos] = spotless_add_sprocedure(prog, vmax_sos, t*(1-t),V_vars,degree-2);
prog = prog.withSOS(vmax_sos);
end

% (3) w(x) >= Wlvl, x in X ==> V >= 0
[prog, v_sos] = spotless_add_sprocedure(prog, V, h_X,V_vars,degree-2);
[prog, v_sos] = spotless_add_sprocedure(prog, v_sos, W-Wlvl,V_vars,degree-W_degree);
[prog, v_sos] = spotless_add_sprocedure(prog, v_sos, t*(1-t),V_vars,degree-2);
prog = prog.withSOS(v_sos);

% % (4) x on bnd(B) ==> V >= 1
% [prog, vbnd_sos] = spotless_add_eq_sprocedure(prog, V-1, B,V_vars,degree-2);
% [prog, vbnd_sos] = spotless_add_sprocedure(prog, vbnd_sos, t*(1-t),V_vars,degree-2);
% prog = prog.withSOS(vbnd_sos);

% % (4) x in X ==> q >= V
% [prog, qbnd_sos] = spotless_add_sprocedure(prog, q-V, h_X,V_vars,degree-2);
% [prog, qbnd_sos] = spotless_add_sprocedure(prog, qbnd_sos, t*(1-t),V_vars,degree-2);
% [prog, qbnd_sos] = spotless_add_sprocedure(prog, qbnd_sos, W,V_vars,degree-W_degree);
% prog = prog.withSOS(qbnd_sos);
% 
% % (4) x in X, W <= 0 ==> q >= 0
% [prog, q_sos] = spotless_add_sprocedure(prog, q, h_X,V_vars,degree-2);
% [prog, q_sos] = spotless_add_sprocedure(prog, q_sos, -W,V_vars,degree-W_degree);
% [prog, q_sos] = spotless_add_sprocedure(prog, q_sos, t*(1-t),V_vars,degree-2);
% prog = prog.withSOS(q_sos);


%% Set up cost function -- integration over a sphere

% cost = 0;
% cost = spotlessIntegral(prog,q,x,R_diag,t,[0 1]);
if mode == 1
  cost = spotlessIntegral(prog,-V,x,R_diag,t,[0 1]);
else
  cost = spotlessIntegral(prog,V,x,R_diag,t,[0 1]);
end
% cost = spotlessIntegral(prog,-q,x,R_diag,[],[]);


spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(cost,solver,spot_options);



%% Plotting
Vsol = subs(sol.eval(V),x,scale.*x);

% qsol = subs(sol.eval(q),x,scale.*x);
% Wsol = subs(sol.eval(W),x,scale.*x);
% model.plotfun(0, sum(Vsol), Vsol, subs(h_X,x,scale.*x), R_diag, msspoly('t',1), x, []);
% sum(Vsol)

if mode == 1
  %% Get level set
  prog = spotsosprog();
  [prog,x]=prog.newIndeterminate('x', model.num_states); % state
  [prog,t]=prog.newIndeterminate('t', 1); % state
  [prog,rho] = prog.newFree(1);
  
  % (2) w(x) = Wlvl, x in X ==> V >= rho
  [prog, vmax_sos] = spotless_add_sprocedure(prog, Vsol - rho, h_X,V_vars,degree-2);
  [prog, vmax_sos] = spotless_add_eq_sprocedure(prog, vmax_sos, W-Wlvl,V_vars,degree-W_degree);
  [prog, vmax_sos] = spotless_add_sprocedure(prog, vmax_sos, t*(1-t),V_vars,degree-2);
  prog = prog.withSOS(vmax_sos);
  
  % (2b) w(x) >= Wlvl, x in dX ==> V >= rho
  [prog, vmax_sos] = spotless_add_eq_sprocedure(prog, Vsol - rho, h_X,V_vars,degree-2);
  [prog, vmax_sos] = spotless_add_sprocedure(prog, vmax_sos, W-Wlvl,V_vars,degree-W_degree);
  [prog, vmax_sos] = spotless_add_sprocedure(prog, vmax_sos, t*(1-t),V_vars,degree-2);
  prog = prog.withSOS(vmax_sos);
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  % solver = @spot_sedumi;
  sol = prog.minimize(-rho,solver,spot_options);
  
  Vlvl = sol.eval(rho);
else
  Vlvl = 1;
end
end