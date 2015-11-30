% Run a 0-step reachability problem (completely passive, no actuation)
checkDependency('spotless')
checkDependency('mosek')

% Setup SOS-program
prog = spotsosprog;

% degree of V,W
degree = 6;

[prog,t]=prog.newIndeterminate('t',1);
[prog,q]=prog.newIndeterminate('q',2);
[prog,v]=prog.newIndeterminate('v',2);

% solve once, then remove cost function and re-solve with cost as
% constraint (to improve numerical conditioning)
do_backoff = false;

% Let V depend on t--probably want it true for this problem class
time_varying = true;


% problem data
T = 2;  %step time. Note that no actual step is taken, just the amount of time allowed to get to the goal region
g = 10; %gravity acceleration
z_nom = 1; %COM height
goal_radius = .01; % Radius of goal region
R_diag = [2 2 2 2];  % state space ball

x = [q;v];

if time_varying
  V_vars = [t;x];
else
  V_vars = x;
end
W_vars = x;

% Create polynomials V(t,x) and W(x)
[prog,V] = prog.newFreePoly(monomials(V_vars,0:degree));
[prog,W] = prog.newFreePoly(monomials(W_vars,0:degree));

% LIPM dynamics
f = [v;q*g/z_nom];

% tau = t/T
% dx/dtau = dx/dt * dt/dtau = dx/dt*T

% time rescaling
f = f*T;
T = 1;


Vdot = diff(V,x)*f + diff(V,t);

% 4 SOS constraints
% (1) V(t,x) >= 0 for x in goal region
% (2) -Vdot(t,x) <= 0 for x in X
% (3) W(x) >= 0 for x in X
% (4) W(x) >= V(0,x) + 1 for x in X
sos = [V;-Vdot; W;W - subs(V,t,0) - 1];

% construct goal region and X ball constraints
A = diag(1./(R_diag.^2));
h_goal = goal_radius^2 - x'*x;
h_X = 1 - x'*A*x;

% add goal region and h_X constraints
[prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), h_goal,V_vars,degree-2);
[prog, sos(2)] = spotless_add_sprocedure(prog, sos(2), h_X,V_vars,degree-2);
[prog, sos(3)] = spotless_add_sprocedure(prog, sos(3), h_X,W_vars,degree-2);
[prog, sos(4)] = spotless_add_sprocedure(prog, sos(4), h_X,W_vars,degree-2);

% 0 <= t < = 1
% could also write this with two constraints
if time_varying
  [prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), T^2-t^2,V_vars,degree-2);
  [prog, sos(2)] = spotless_add_sprocedure(prog, sos(2), T^2-t^2,V_vars,degree-2);
end

%% Setup cost function
[vars,alphas,coeff] = decomp(W, prog.freeVar);
inds = [1;3;2;4];
alphas = alphas(:,inds);
assert(isequal(vars(inds),x))

nX = 2;

% [~,alphas] = monomials(x_vars,0:degree);
betas = 0.5*(alphas + 1);
Ra = (1.^(sum(alphas,2) + nX))./(sum(alphas,2) + nX);
IS = 2*prod(gamma(betas),2)./(gamma(sum(betas,2)));
l = Ra.*IS;
alphaszero = (mod(alphas,2) ~= 0);
alphaszero = any(alphaszero,2);
l(alphaszero) = 0;
l = l.*prod(repmat(R_diag,size(alphas,1),1).^(alphas+1),2);

for i=1:length(sos)
  prog = prog.withSOS(sos(i));
end

options = spotprog.defaultOptions;
options.verbose = true;
options.do_fr = true;
sol = prog.minimize(coeff*l,@spot_mosek,options);

if do_backoff
  % resolve problem with cost replaced by a constraint
  prog = prog.withPos(sol.eval(coeff*l)*1.01 - coeff*l);
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

figure(2)
hold off
contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0],{'b','r'});
xlabel('q_1')
ylabel('v_1')
title('V(0,x)')

%%
V0 = Vsol;
save V0_LIPM V0
