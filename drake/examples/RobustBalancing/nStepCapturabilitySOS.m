function nStepCapturabilitySOS(model, T, R_diag, target, n, options)
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
if ~isfield(options,'scale_input')
  scale_input = ones(model.num_inputs,1);
elseif length(options.scale_input) == 1
  scale_input = ones(model.num_inputs,1)*options.scale_input;
else
  scale_input = options.scale_input;
end
scale_input_inv = 1./scale_input;

%% Solution method settings
degree = options.degree; % degree of V,W
time_varying = n > 0 || model.num_inputs; % Let V depend on t--probably want it true for this problem class

%% Create SOS program
prog = spotsosprog;

%% Create indeterminates
[prog,t]=prog.newIndeterminate('t',1); % time
[prog,x]=prog.newIndeterminate('x', model.num_states); % state
if model.num_inputs > 0
  [prog,u]=prog.newIndeterminate('u', model.num_inputs); % input
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
  filename = solutionFileName(model, n - 1);
  if ~exist(filename, 'file')
    display(sprintf('Solving for %d-step first',n-1))
    nStepCapturabilitySOS(model, T, R_diag, target, n - 1, options);
  end
  data = load(filename);
  V0 = subs(data.Vsol,x,x.*scale_inv);
end

%% Scale r_diag
R_diag = scale'.*R_diag;

%% Create polynomials V(t,x) and W(x)
if time_varying
  V_vars = [t;x];
else
  V_vars = x;
end
W_vars = x;
[prog,V] = prog.newFreePoly(monomials(V_vars,0:degree));
[prog,W] = prog.newFreePoly(monomials(W_vars,0:degree));

%% Dynamics
f = scale.*model.dynamics(t, scale_inv.*x, scale_input_inv.*u);

% Time rescaling
% tau = t / T
% dx/dtau = dx/dt * dt/dtau = dx/dt*T
T_init = T;
f = f*T;
T = 1;

Vdot = diff(V,x)*f + diff(V,t);

Vdot_degree = even_degree(Vdot,[x;u]);


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
sos = msspoly;
if n > 0
  % (1) V(T,x) >= 0 for x in goal region
  % goal region
  if options.free_final_time
    V_goal_eqn = V*(1+[V_vars;s]'*[V_vars;s]);
    goal_vars = [V_vars;s];
  else
    V_goal_eqn = subs(V,t,T)*(1+[x;s]'*[x;s]);
    goal_vars = [W_vars;s];
  end
  [prog, goal_sos] = spotless_add_sprocedure(prog, V_goal_eqn, V0p,goal_vars,2);

  % state constraint
  [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, h_X,goal_vars,degree);

  % reset map input limits
  [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, model.resetInputLimits(s),goal_vars,degree);
  
  if options.free_final_time && time_varying
    [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, t * (T - t),goal_vars,degree);
  end
  
  sos = [sos; goal_sos];
else
  if ~isempty(target)
    % (1) V(t,x) >= 0 for x in goal region
    [prog, goal_sos] = spotless_add_sprocedure(prog, V, V0p,V_vars,degree-2);
    
    if time_varying
      [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, t * (T - t),V_vars,degree-2);
    end
    
    sos = [sos; goal_sos];
  else
    if options.free_final_time
      [prog, goal_sos] = spotless_add_sprocedure(prog, subs(V,x,zeros(model.num_states,1)), t*(T-t),t,degree-2);
      sos = [sos; goal_sos];
    else
      prog = prog.withPos(subs(subs(V,t,T),x,zeros(model.num_states,1)));      
    end
  end
end


% (2) -Vdot(t,x,u) <= 0 for x in X
[prog, Vdot_sos] = spotless_add_sprocedure(prog, -Vdot, h_X,[V_vars;u],Vdot_degree-2);

% input limits
[prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, model.inputLimits(scale_input_inv.*u, scale_inv.*x),[V_vars;u],[]);

input_equality_constraints = model.inputEqualityConstraints(scale_input_inv.*u, scale_inv.*x);
input_equality_constraint_degree = even_degree(input_equality_constraints,[x;u]);
for i = 1 : length(input_equality_constraints) % TODO iteration in spotless_add_eq_sprocedure
  [prog, Vdot_sos] = spotless_add_eq_sprocedure(prog, Vdot_sos, input_equality_constraints(i), [V_vars; u], input_equality_constraint_degree); % TODO: degree
end

% 0 <= t < = T
% could also write this with two constraints
if time_varying
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, t * (T - t),[V_vars;u],Vdot_degree-2);
end
sos = [sos; Vdot_sos];

% (3) W(x) >= 0 for x in X
[prog, sos(end + 1)] = spotless_add_sprocedure(prog, W, h_X,W_vars,degree-2);

% (4) W(x) >= V(0,x) + 1 for x in X
[prog, sos(end + 1)] = spotless_add_sprocedure(prog, W - subs(V,t,0) - 1, h_X,W_vars,degree-2);

for i=1:length(sos)
  prog = prog.withSOS(sos(i));
end

%% Set up cost function -- integration over a sphere
cost = spotlessIntegral(prog,W,x,R_diag,[],[]);

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

%% Plotting
Vsol = subs(sol.eval(V),x,scale.*x);
Wsol = subs(sol.eval(W),x,scale.*x);
R_diag = scale_inv'.*R_diag;
model.plotfun(n, Vsol, Wsol, subs(h_X,x,scale.*x), R_diag, t, x);

%%
T = T_init;
save(solutionFileName(model, n),'Vsol','model','T','R_diag')

end

function filename = solutionFileName(model, n)
filename_suffix = class(model);
filename = sprintf(['V%d_' filename_suffix '.mat'], n);
end
