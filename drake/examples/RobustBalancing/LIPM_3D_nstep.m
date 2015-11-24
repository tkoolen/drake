clear all

n = 2; % step count

data=load(sprintf('V%d_LIPM',n-1));
V0 = data.Vsol;
prog = spotsosprog;
degree = 6;

do_backoff = false;
time_varying = true;

T = .3;
g = 10;
z_nom = 1;
step_max = .7;

[prog,t]=prog.newIndeterminate('t',1);
[prog,q]=prog.newIndeterminate('q',2);
[prog,v]=prog.newIndeterminate('v',2);
[prog,u]=prog.newIndeterminate('u',2);
x = [q;v];


if time_varying
  V_vars = [t;x];
else
  V_vars = x;
end

W_vars = x;

[prog,V] = prog.newFreePoly(monomials(V_vars,0:degree));
[prog,W] = prog.newFreePoly(monomials(W_vars,0:degree));

f = [v;q*g/z_nom];

% time rescaling
T_unscaled = T;
f = f*T;
T = 1;

Vdot = diff(V,x)*f + diff(V,t);


% from Koolen et. al IJRR
% regions should depend on the instantaneous capture point
r_ic = q + v*sqrt(z_nom/g);

%% jump equation
% qp = -[qm(2);qm(1)]
% vp = [vm(1);-vm(1)];
xp = [q-u;v];
V0p = subs(V0,[x;t],[xp;0]);

sos = [(subs(V,t,T))*(1+x'*x + u'*u);-Vdot; W;W - subs(V,t,0) - 1];

R_diag = [2 2 2 2];
A = diag(1./(R_diag.^2));

h_X = 1 - x'*A*x;
% h_X0 = 1 - q(1)^2 - v(1)^2;




%%Constraints
[prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), V0p,[W_vars;u],2);
[prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), h_X,[W_vars;u],degree);
% [prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), subs(h_X0,x,xp),V_vars,degree);

[prog, sos(2)] = spotless_add_sprocedure(prog, sos(2), h_X,V_vars,degree-2);
[prog, sos(3)] = spotless_add_sprocedure(prog, sos(3), h_X,W_vars,degree-2);
[prog, sos(4)] = spotless_add_sprocedure(prog, sos(4), h_X,W_vars,degree-2);

[prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), step_max^2-u'*u,[W_vars;u],degree);

if time_varying
  [prog, sos(2)] = spotless_add_sprocedure(prog, sos(2), T^2-t^2,V_vars,degree-2);
end
%% Setup cost function
[vars,alphas,coeff] = decomp(W, prog.freeVar);
inds = [1;3;2;4];
alphas = alphas(:,inds);
assert(isequal(vars(inds),x))

nX = length(x);

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


% prog = prog.withPos(850 - coeff*l);
% l = l*0;

options = spotprog.defaultOptions;
options.verbose = true;
options.do_fr = true;
sol = prog.minimize(coeff*l,@spot_mosek,options);

if do_backoff
  
  prog = prog.withPos(sol.eval(coeff*l)*1.01 - coeff*l);
  sol = prog.minimize(0,@spot_mosek,options);
end

%% New-style plotting
Vsol = sol.eval(V);
Wsol = sol.eval(W);
sub_vars = [q(2);v(2);t];
sub_val = [0;0;0];
plot_vars = [q(1);v(1)];

figure(1)
contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[1 0],{'b','r'});

figure(2)
hold off
contourSpotless([Vsol;h_X;r_ic'*r_ic],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0 .07],{'b','r','g'});
%%
save V1_LIPM Vsol