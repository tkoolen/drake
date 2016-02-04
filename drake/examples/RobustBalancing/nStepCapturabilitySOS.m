function [Vsol,Wsol,u_sol] = nStepCapturabilitySOS(model, T, R_diag, target, n, options)
% Run an n-step reachability problem
% @param R_diag state space ball

% TODO: put more stuff in options struct

checkDependency('spotless');
checkDependency('mosek');

if ~isfield(options,'do_backoff')
  options.do_backoff = false;
end

if ~isfield(options,'backoff_ratio')
  options.backoff_ratio = 1.01;
end

if ~isfield(options,'free_final_time')
  options.free_final_time = false; % goal region is not restricted to t=T if true
end

if ~isfield(options,'control_design')
  options.control_design = false;
end

if ~isfield(options,'korda_control_design')
  options.korda_control_design = false;
end

if ~isfield(options,'beta')
  options.beta = 0;
end

% then don't use W at all
if ~isfield(options,'infinite_time')
  options.infinite_time = false;
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

%scaling of input vector
if ~isfield(options,'scale_input') || options.control_design
  scale_input = ones(model.num_inputs,1);
elseif length(options.scale_input) == 1
  scale_input = ones(model.num_inputs,1)*options.scale_input;
else
  scale_input = options.scale_input;
end
scale_input_inv = 1./scale_input;

%% Solution method settings
degree = options.degree; % degree of V,W
time_varying = (n > 0 || model.num_inputs) && ~options.infinite_time; % Let V depend on t--probably want it true for this problem class

%% Create SOS program
prog = spotsosprog;

%% Create indeterminates
if time_varying
  [prog,t]=prog.newIndeterminate('t',1); % time
else
  t = msspoly('t',1);
end
[prog,x]=prog.newIndeterminate('x', model.num_states); % state
if model.num_inputs > 0
  if ~options.control_design
    [prog,u]=prog.newIndeterminate('u', model.num_inputs); % input
  else
    u = msspoly('u',model.num_inputs);
  end
else
%   u = msspoly;
  u = zeros(0,1);
end
if n > 0
  if model.num_reset_inputs > 0
    [prog,s]=prog.newIndeterminate('s', model.num_reset_inputs); % reset map input
  else
    s = [];
  end
end

%% Load previous problem data
if n > 0
  if isfield(options,'V0')
    V0 = subs(options.V0,x,x.*scale_inv);
  else
    filename = solutionFileName(model, n - 1);
    if ~exist(filename, 'file')
      display(sprintf('Solving for %d-step first',n-1))
      nStepCapturabilitySOS(model, T, R_diag, target, n - 1, options);
    end
    data = load(filename);
    V0 = subs(data.Vsol,x,x.*scale_inv);
  end
end

%% Scale r_diag
R_diag = scale'.*R_diag;

%% Create polynomials V(t,x) and W(x)
if time_varying
  V_vars = [t;x];
else
  V_vars = x;
end
[prog,V] = prog.newFreePoly(monomials(V_vars,0:degree));

if ~options.infinite_time
  W_vars = x;
  [prog,W] = prog.newFreePoly(monomials(W_vars,0:degree));
else
  W_vars = V_vars;
end

%% Dynamics
f = scale.*model.dynamics(t, scale_inv.*x, scale_input_inv.*u);

% Time rescaling
% tau = t / T
% dx/dtau = dx/dt * dt/dtau = dx/dt*T

T_init = T;
f = f*T;
T = 1;

Vdot = diff(V,x)*f + diff(V,t);

if options.control_design && model.num_inputs > 0
  Vdot_degree = even_degree(Vdot,x);
  
  [prog,p] = prog.newFreePoly(monomials(V_vars,0:Vdot_degree),model.num_inputs);
  [A_u,b_u,C_u,d_u] = model.unitBoxInputTransform();
  Vdot = subs(Vdot,u,C_u*u + d_u);
  
  dVdotdu = diff(Vdot,u);
  if deg(Vdot,u) ~= 1
    error('System must be control affine');
  end
  
  if options.korda_control_design
    Vdot = subs(Vdot,u,-ones(model.num_inputs,1)) + sum(p) - options.beta*V; % u <--- (u + 1)/2, plus setting u_bar = 1 (from paper).
    dVdotdu = 2*dVdotdu;
  else
    Vdot = subs(Vdot,u,u*0) + sum(p);
  end
else
  Vdot_degree = even_degree(Vdot,[x;u]);
end


%% Goal region
if n > 0
  % jump equation
  xp = scale.*model.reset(t, scale_inv.*x, s);

  % for n > 0, goal region is based off V from 0-step model
  % V0p(x) = V0(0,xp)
  V0p = subs(V0,[x;t],[xp;0]);
else
  if ~isempty(target)
    V0p = target(scale_inv.*x);
  end
end

% State constraint
A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x;

%% SOS constraints
if options.infinite_time
  V_goal_min = 1;
else
  V_goal_min = 0;
end

if n > 0
  % (1) V(T,x) >= V_goal_min for x in goal region
  % goal region
  if options.free_final_time
    V_goal_eqn = (V-V_goal_min)*(1+[V_vars;s]'*[V_vars;s]);
    goal_vars = [V_vars;s];
  else
    V_goal_eqn = (subs(V,t,T)-V_goal_min)*(1+[x;s]'*[x;s]);
    goal_vars = [W_vars;s];
  end
  [prog, goal_sos] = spotless_add_sprocedure(prog, V_goal_eqn, V0p,goal_vars,degree);

  % state constraint
  [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, h_X,goal_vars,degree);

  % reset map input limits
  [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, model.resetInputLimits(s),goal_vars,degree);
  
  if options.free_final_time && time_varying
    [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, t * (T - t),goal_vars,degree);
  end
  
  prog = prog.withSOS(goal_sos);
else
  if ~isempty(target)
    % (1) V(t,x) >= 0 for x in goal region
    [prog, goal_sos] = spotless_add_sprocedure(prog, V-V_goal_min, V0p,V_vars,degree-2);
    
    if time_varying
      [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, t * (T - t),V_vars,degree-2);
    end
    
    prog = prog.withSOS(goal_sos);
  else
    if options.free_final_time
      [prog, goal_sos] = spotless_add_sprocedure(prog, subs(V-V_goal_min,x,zeros(model.num_states,1)), t*(T-t),t,degree-2);
      prog = prog.withSOS(goal_sos);
    else
      prog = prog.withPos(subs(subs(V-V_goal_min,t,T),x,zeros(model.num_states,1)));
    end
  end
end


% (2) -Vdot(t,x,u) <= 0 for x in X
if options.control_design
  Vdot_vars = V_vars;
else
  Vdot_vars = [V_vars;u];
end

[prog, Vdot_sos] = spotless_add_sprocedure(prog, -Vdot, h_X,Vdot_vars,Vdot_degree-2);

if ~options.control_design
  % input limits
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, model.inputLimits(scale_input_inv.*u, scale_inv.*x),Vdot_vars,[]);
  
  input_equality_constraints = model.inputEqualityConstraints(scale_input_inv.*u, scale_inv.*x);
  input_equality_constraint_degree = even_degree(input_equality_constraints,[x;u]);
  
  for i = 1 : length(input_equality_constraints) % TODO iteration in spotless_add_eq_sprocedure
    [prog, Vdot_sos] = spotless_add_eq_sprocedure(prog, Vdot_sos, input_equality_constraints(i), Vdot_vars, input_equality_constraint_degree); % TODO: degree
  end
end

% 0 <= t < = T
% could also write this with two constraints
if time_varying
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t * (T - t),Vdot_vars,Vdot_degree-2);
end
[prog,Vdot_ind] = prog.withSOS(Vdot_sos);

if options.control_design
    p_pos_sos = msspoly;
    p_neg_sos = msspoly;
    p_pos_ind = [];
    p_neg_ind = [];
    p_sos_ind = [];
    for i=1:model.num_inputs,
      p_pos_sos_i = p(i) - dVdotdu(i);
      [prog, p_pos_sos_i] = spotless_add_sprocedure(prog, p_pos_sos_i, h_X,V_vars,Vdot_degree-2);
      if time_varying
        [prog, p_pos_sos_i] = spotless_add_sprocedure(prog, p_pos_sos_i,  t * (T - t),V_vars,Vdot_degree-2);
      end
      [prog,p_pos_ind(i)] = prog.withSOS(p_pos_sos_i);
      p_pos_sos = [p_pos_sos;p_pos_sos_i];

      if options.korda_control_design
        p_neg_sos_i = p(i);
      else
        p_neg_sos_i = p(i) + dVdotdu(i);
      end
      [prog, p_neg_sos_i] = spotless_add_sprocedure(prog, p_neg_sos_i, h_X,V_vars,Vdot_degree-2);
      if time_varying
        [prog, p_neg_sos_i] = spotless_add_sprocedure(prog, p_neg_sos_i,  t * (T - t),V_vars,Vdot_degree-2);
      end
      [prog,p_neg_ind(i)] = prog.withSOS(p_neg_sos_i);

      p_neg_sos = [p_neg_sos;p_neg_sos_i];

      if ~options.korda_control_design
        [prog, p_sos_i] = spotless_add_sprocedure(prog, p, h_X,V_vars,Vdot_degree-2);
        if time_varying
          [prog, p_sos_i] = spotless_add_sprocedure(prog, p_sos_i,  t * (T - t),V_vars,Vdot_degree-2);
        end
        [prog,p_sos_ind(i)] = prog.withSOS(p_sos_i);
      end
    end
end
  

if options.infinite_time
  % (3) V(x) >= -1 for x in X
  [prog, V_min_sos] = spotless_add_sprocedure(prog, V+1, h_X,V_vars,degree-2);
  prog = prog.withSOS(V_min_sos);
else
  % (3) W(x) >= 0 for x in X
  [prog, W_sos] = spotless_add_sprocedure(prog, W, h_X,W_vars,degree-2);
  [prog, W_ind] = prog.withSOS(W_sos);

  % (4) W(x) >= V(0,x) + 1 for x in X
  [prog, WminusV_sos] = spotless_add_sprocedure(prog, W - subs(V,t,0) - 1, h_X,W_vars,degree-2);
  [prog, WminusV_ind] = prog.withSOS(WminusV_sos);
end
%% Set up cost function -- integration over a sphere

if options.infinite_time
  cost = spotlessIntegral(prog,V,x,R_diag,[],[]);
else
  cost = spotlessIntegral(prog,W,x,R_diag,[],[]);
end

%% Solve
spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(cost,solver,spot_options);

if options.do_backoff
  % resolve problem with cost replaced by a constraint
  prog = prog.withPos(sol.eval(cost)*options.backoff_ratio - cost);
  sol = prog.minimize(0,solver,spot_options);
end

%% controller extraction (Majumdar et al. Mark's old code)
if options.control_design
  % get sos decomps  
  y_vdot = sol.prog.sosEqsDualVars{Vdot_ind};
  basis_vdot = sol.prog.sosEqsBasis{Vdot_ind};
  mu_vdot = double(sol.dualEval(y_vdot));
  [M_vdot,G] = momentMatrix(mu_vdot,basis_vdot);
  
  M_sig = cell(model.num_inputs,1);
  u_sol = msspoly;
  for i=1:model.num_inputs,
    y_p_pos = sol.prog.sosEqsDualVars{p_pos_ind(i)};
    basis_p_pos = sol.prog.sosEqsBasis{p_pos_ind(i)};

    y_p_neg = sol.prog.sosEqsDualVars{p_neg_ind(i)};
    basis_p_neg = sol.prog.sosEqsBasis{p_neg_ind(i)};
    
    if length(basis_vdot) ~= length(basis_p_pos)
      err1 = 1;
    elseif length(basis_p_neg) ~= length(basis_p_pos)
      err2 = 1;
    else
      [~,err1] = double(basis_vdot-basis_p_pos);
      [~,err2] = double(basis_p_neg-basis_p_pos);
    end
    if err1 || err2
      warning(['Same basis did not come out of SOS ' ...
        'Decomposition.']);
      keyboard
    end
    
    if options.korda_control_design
      y_sig = double(sol.dualEval(y_p_pos));
    else
      y_sig = double(sol.dualEval(y_p_pos-y_p_neg));
    end
    M_sig{i} = momentMatrix(y_sig,basis_vdot);
    
    [u_i,u_i_coeff] = solveController(M_vdot,M_sig{i},G);
    u_sol = [u_sol;u_i];
    uMomentBasis{i} = G;
    uCoeff{i} = u_i_coeff;
  end
  u_sol = C_u*u_sol + d_u;
  u_sol = subs(u_sol,x,scale.*x);

  if options.korda_control_design
    u_sol = 2*u_sol - 1;
  end
  clean(u_sol,1e-4)
end

%% Plotting
Vsol = subs(sol.eval(V),x,scale.*x);
if options.infinite_time
  Wsol = subs(sol.eval(V),x,scale.*x);
else
  Wsol = subs(sol.eval(W),x,scale.*x);
end
R_diag = scale_inv'.*R_diag;
if options.control_design
  model.plotfun(n, Vsol, Wsol, subs(h_X,x,scale.*x), R_diag, t, x, u_sol);
else
  model.plotfun(n, Vsol, Wsol, subs(h_X,x,scale.*x), R_diag, t, x, []);
end

%%
T = T_init;
Vsol = subs(Vsol,t,t/T);
if options.control_design
  u_sol = subs(u_sol,t,t/T);
end

vars_to_save = {'t', 'x', 'Vsol', 'model', 'T', 'R_diag'};
if options.control_design
  save(solutionFileName(model, n),'Vsol','model','T','R_diag','u_sol','uMomentBasis', 'uCoeff')
else
  save(solutionFileName(model, n),'Vsol','model','T','R_diag')
end
save(solutionFileName(model, n), vars_to_save{:});
end

function filename = solutionFileName(model, n)
filename_suffix = class(model);
filename = sprintf(['V%d_' filename_suffix '.mat'], n);
end
