function LIPM_3D_nstep(n)
% Run an n-step reachability problem for the LIPM
% Constant height, angular momentum model
% Control input is the foot position on each step (massless foot)

checkDependency('spotless');
checkDependency('mosek');


%% Load previous problem data
if n > 0
  filename = sprintf('V%d_LIPM.mat',n-1);
  if ~exist(filename, 'file')
    LIPM_3D_nstep(n - 1);
  end
  data=load(sprintf('V%d_LIPM',n-1));
  V0 = data.Vsol;
end

%% Solution method settings
degree = 6; % degree of V,W
do_backoff = false; % solve once, then remove cost function and re-solve with cost as constraint (to improve numerical conditioning)
time_varying = true; % Let V depend on t--probably want it true for this problem class
R_diag = [2 2 2 2]; % state space ball

%% Model parameters
if n > 0
  T = 0.3;  % Step-time
else
  T = 2;
end
g = 10;  % gravity acceleration
z_nom = 1; % nominal center of mass height
step_max = .7; % max step distance

%% Create SOS program
prog = spotsosprog;

%% Create indeterminates
[prog,t]=prog.newIndeterminate('t',1);
[prog,q]=prog.newIndeterminate('q',2);
[prog,v]=prog.newIndeterminate('v',2);
if n > 0
  [prog,u]=prog.newIndeterminate('u',2);
end
x = [q;v];

%% Create polynomials V(t,x) and W(x)
if time_varying
  V_vars = [t;x];
else
  V_vars = x;
end
W_vars = x;
[prog,V] = prog.newFreePoly(monomials(V_vars,0:degree));
[prog,W] = prog.newFreePoly(monomials(W_vars,0:degree));

%% LIPM dynamics
f = [v;q*g/z_nom];

% Time rescaling
% tau = t / T
% dx/dtau = dx/dt * dt/dtau = dx/dt*T
T_unscaled = T;
f = f*T;
T = 1;

Vdot = diff(V,x)*f + diff(V,t);

%% Goal region
if n > 0
  % jump equation
  % control input changes q only
  % qp = qm - u
  xp = [q-u;v];

  % for n > 0, goal region is based off V from 0-step model
  % V0p(x) = V0(0,xp)
  V0p = subs(V0,[x;t],[xp;0]);
else
  % use a small radius around the origin
  goal_radius = .01; % TODO: setting
  V0p = goal_radius^2 - x'*x;
end

% State constraint
A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x;

%% SOS constraints
sos = msspoly;
if n > 0
  % (1) V(T,x) >= 0 for x in goal region
  % goal region
  [prog, goal_sos] = spotless_add_sprocedure(prog, (subs(V,t,T))*(1+x'*x + u'*u), V0p,[W_vars;u],2);

  % state constraint
  [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, h_X,[W_vars;u],degree);

  % control input limit -step_max <= u <= step_max
  [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, step_max^2-u'*u,[W_vars;u],degree);
else
  % (1) V(t,x) >= 0 for x in goal region
  [prog, goal_sos] = spotless_add_sprocedure(prog, V, V0p,V_vars,degree-2);

  if time_varying
    [prog, goal_sos] = spotless_add_sprocedure(prog, goal_sos, T^2-t^2,V_vars,degree-2);
  end
end
sos = [sos; goal_sos];


% (2) -Vdot(t,x) <= 0 for x in X
[prog, Vdot_sos] = spotless_add_sprocedure(prog, -Vdot, h_X,V_vars,degree-2);
% 0 <= t < = T
% could also write this with two constraints
if time_varying
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, T^2-t^2,V_vars,degree-2);
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
options = spotprog.defaultOptions;
options.verbose = true;
options.do_fr = true;
sol = prog.minimize(cost,@spot_mosek,options);

if do_backoff
  % resolve problem with cost replaced by a constraint
  prog = prog.withPos(sol.eval(cost)*1.01 - cost);
  sol = prog.minimize(0,@spot_mosek,options);
end

%% Plotting
Vsol = sol.eval(V);
Wsol = sol.eval(W);
sub_vars = [q(2);v(2);t];
sub_val = [0;0;0];
plot_vars = [q(1);v(1)];

figure(1)
contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[1 0],{'b','r'});
xlabel('q_1')
ylabel('v_1')
title('W(x)')


% from Koolen et. al IJRR
% regions should depend on the instantaneous capture point
r_ic = q + v*sqrt(z_nom/g);
dN = captureLimit(T_unscaled, 0, step_max, z_nom, g, n); % theoretical max ICP distance

figure(n*10+2)
contourSpotless([Vsol;h_X;r_ic'*r_ic],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0 dN^2],{'b','r','g'});
xlabel('q_1')
ylabel('v_1')
title('V(0,x)')

%%
save(sprintf('V%d_LIPM',n),'Vsol')

end
