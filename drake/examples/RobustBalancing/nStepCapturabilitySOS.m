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

%% Solution method settings
degree = options.degree; % degree of V,W
time_varying = n > 0 || model.num_inputs; % Let V depend on t--probably want it true for this problem class

%% Load previous problem data
if n > 0
  filename = solutionFileName(model, n - 1);
  if ~exist(filename, 'file')
    nStepCapturabilitySOS(model, T, R_diag, target, n - 1, options);
  end
  data = load(filename);
  V0 = data.Vsol;
end

%% Create SOS program
prog = spotsosprog;

%% Create indeterminates
[prog,t]=prog.newIndeterminate('t',1); % time
[prog,x]=prog.newIndeterminate('x', model.num_states); % state
if model.num_inputs > 0
  [prog,u]=prog.newIndeterminate('u', model.num_inputs); % input
else
  u = msspoly;
end
if n > 0
  [prog,s]=prog.newIndeterminate('s', model.num_reset_inputs); % reset map input
end

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
f = model.dynamics(t, x, u);

% Time rescaling
% tau = t / T
% dx/dtau = dx/dt * dt/dtau = dx/dt*T
f = f*T;
T = 1;

Vdot = diff(V,x)*f + diff(V,t);

Vdot_degree = even_degree(Vdot,[x;u]);


%% Goal region
if n > 0
  % jump equation
  xp = model.reset(t, x, s);

  % for n > 0, goal region is based off V from 0-step model
  % V0p(x) = V0(0,xp)
  V0p = subs(V0,[x;t],[xp;0]);
else
  if ~isempty(target)
    V0p = target(x);
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
  [prog, goal_sos] = spotless_add_sprocedure(prog, (subs(V,t,T))*(1+x'*x + s'*s), V0p,[W_vars;s],2);

  % state constraint
  [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, h_X,[W_vars;s],degree);

  % reset map input limits
  [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, model.resetInputLimits(s),[W_vars;s],degree);
  
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
    prog = prog.withPos(subs(subs(V,t,T),x,zeros(model.num_states,1)));
  end
end


% (2) -Vdot(t,x,u) <= 0 for x in X
[prog, Vdot_sos] = spotless_add_sprocedure(prog, -Vdot, h_X,[V_vars;u],Vdot_degree-2);

% input limits
input_limit_degree = even_degree(model.inputLimits(u,x),[x;u]);
[prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, model.inputLimits(u, x),[V_vars;u],[]);

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
sol = prog.minimize(cost,@spot_mosek,spot_options);

if options.do_backoff
  % resolve problem with cost replaced by a constraint
  prog = prog.withPos(sol.eval(cost)*options.backoff_ratio - cost);
  sol = prog.minimize(0,@spot_mosek,spot_options);
 end

%% Plotting
Vsol = sol.eval(V);
Wsol = sol.eval(W);
model.plotfun(n, Vsol, Wsol, h_X, R_diag, t, x);

%%
save(solutionFileName(model, n),'Vsol')

end

function filename = solutionFileName(model, n)
filename_suffix = class(model);
filename = sprintf(['V%d_' filename_suffix '.mat'], n);
end
