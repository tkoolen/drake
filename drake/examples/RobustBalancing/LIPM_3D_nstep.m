function LIPM_3D_nstep(n)
% Run an n-step reachability problem for the LIPM
% Constant height, angular momentum model
% Control input is the foot position on each step (massless foot)
checkDependency('spotless');
checkDependency('mosek');

if n == 0
  LIPM_3D_zerostep();
  return;
end

% load previous problem data
filename = sprintf('V%d_LIPM.mat',n-1);
if ~exist(filename, 'file')
  LIPM_3D_nstep(n - 1);
end
data=load(sprintf('V%d_LIPM',n-1));
V0 = data.Vsol;

% solution method settings
degree = 6; % degree of V,W
do_backoff = false; % solve once, then remove cost function and re-solve with cost as constraint (to improve numerical conditioning)
time_varying = true; % Let V depend on t--probably want it true for this problem class
R_diag = [2 2 2 2]; % state space ball

% model parameters
T = .3;  % Step-time
g = 10;  % gravity acceleration
z_nom = 1; % nominal center of mass height
step_max = .7; % max step distance

% Create SOS program
prog = spotsosprog;

[prog,t]=prog.newIndeterminate('t',1);
[prog,q]=prog.newIndeterminate('q',2);
[prog,v]=prog.newIndeterminate('v',2);
[prog,u]=prog.newIndeterminate('u',2);
x = [q;v];


% Create polynomials V(t,x) and W(x)
if time_varying
  V_vars = [t;x];
else
  V_vars = x;
end
W_vars = x;
[prog,V] = prog.newFreePoly(monomials(V_vars,0:degree));
[prog,W] = prog.newFreePoly(monomials(W_vars,0:degree));

% LIPM dynamics
f = [v;q*g/z_nom];

% time rescaling
% dx/dtau = dx/dt * dt/dtau = dx/dt*T
T_unscaled = T;
f = f*T;
T = 1;

Vdot = diff(V,x)*f + diff(V,t);

% goal region
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


A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x;

% 4 SOS constraints
% (1) V(T,x) >= 0 for x in goal region
% (2) -Vdot(t,x) <= 0 for x in X
% (3) W(x) >= 0 for x in X
% (4) W(x) >= V(0,x) + 1 for x in X
sos = [(subs(V,t,T))*(1+x'*x + u'*u);-Vdot; W;W - subs(V,t,0) - 1];

% add goal region and h_X constraints
[prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), V0p,[W_vars;u],2);
[prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), h_X,[W_vars;u],degree);
[prog, sos(2)] = spotless_add_sprocedure(prog, sos(2), h_X,V_vars,degree-2);
[prog, sos(3)] = spotless_add_sprocedure(prog, sos(3), h_X,W_vars,degree-2);
[prog, sos(4)] = spotless_add_sprocedure(prog, sos(4), h_X,W_vars,degree-2);

% control input limit -step_max <= u <= step_max
[prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), step_max^2-u'*u,[W_vars;u],degree);

% 0 <= t < = 1
% could also write this with two constraints
if time_varying
  [prog, sos(2)] = spotless_add_sprocedure(prog, sos(2), T^2-t^2,V_vars,degree-2);
end
%% Setup cost function -- integration over a sphere and then solve
cost = spotlessIntegral(prog,W,x,R_diag,[],[]);

for i=1:length(sos)
  prog = prog.withSOS(sos(i));
end

options = spotprog.defaultOptions;
options.verbose = true;
options.do_fr = true;
sol = prog.minimize(cost,@spot_mosek,options);

if do_backoff
  % resolve problem with cost replaced by a constraint
  prog = prog.withPos(sol.eval(cost)*1.01 - cost);
  sol = prog.minimize(0,@spot_mosek,options);
end

%% plotting
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
