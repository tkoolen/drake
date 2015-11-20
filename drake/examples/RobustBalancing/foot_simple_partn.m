clear all

ind = 2;
filename = sprintf('V%d',ind-1);
Vprev = load(filename);
V0 = Vprev.Vsol;

prog = spotsosprog;
degree = 6;

do_backoff = false;
time_varying = true;

T = 1;
g = 10;
z_nom = 1;

[prog,t]=prog.newIndeterminate('t',1);
[prog,q]=prog.newIndeterminate('q',2);
[prog,v]=prog.newIndeterminate('v',2);
[prog,u]=prog.newIndeterminate('u',1);
x = [q;v];


if time_varying
  V_vars = [t;x];
else
  V_vars = x;
end

W_vars = x;

[prog,V] = prog.newFreePoly(monomials(V_vars,0:degree));
[prog,W] = prog.newFreePoly(monomials(W_vars,0:degree));


f = [v;q(1)*g/z_nom;u];

% time rescaling
T_unscaled = T;
f = f*T;
T = 1;

Vdot = diff(V,x)*f + diff(V,t);
%% jump equation
% qp = -[qm(2);qm(1)]
% vp = [vm(1);-vm(1)];
xp = [-q(2);-q(1);v(1);-v(1)];
V0p = subs(V0,[x;t],[xp;0]);

sos = [(V)*(1+x'*x);-Vdot; W;W - subs(V,t,0) - 1];

R_diag = [2 2 2 2];
A = diag(1./(R_diag.^2));

h_X = 1 - x'*A*x;




%%Constraints
[prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), V0p,x,2);
[prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), h_X,x,degree);
[prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), subs(h_X,x,xp),x,degree);

[prog, sos(2)] = spotless_add_sprocedure(prog, sos(2), h_X,[x;u],degree-2);
[prog, sos(3)] = spotless_add_sprocedure(prog, sos(3), h_X,x,degree-2);
[prog, sos(4)] = spotless_add_sprocedure(prog, sos(4), h_X,x,degree-2);

[prog, sos(2)] = spotless_add_sprocedure(prog, sos(2), 1-u^2,[x;u],degree-2);

if time_varying
  [prog, sos(1)] = spotless_add_sprocedure(prog, sos(1), T^2-t^2,V_vars,degree);
  [prog, sos(2)] = spotless_add_sprocedure(prog, sos(2), T^2-t^2,V_vars,degree-2);
end

%% Setup cost function
[vars,alphas,coeff] = decomp(W, prog.freeVar);
vars = vars([1 3 2 4]);
alphas = alphas(:,[1 3 2 4]);
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

% prog = prog.withPos(40 - coeff*l);
% l = l*0;

options = spotprog.defaultOptions;
options.verbose = 1;
options.do_fr = true;
sol = prog.minimize(coeff*l,@spot_mosek,options);

if do_backoff
  
  prog = prog.withPos(sol.eval(coeff*l)*1.01 - coeff*l);
  sol = prog.minimize(0,@spot_mosek,options);
end

%% plotting
Vsol = sol.eval(V);
Wsol = sol.eval(W);

sub_vars = [q(2);v(2);t];
sub_val = [0;0;0];

[Q,QD] = meshgrid(linspace(-R_diag(1),R_diag(1),400),linspace(-R_diag(3),R_diag(3),400));
h_goal_val = reshape(msubs(subs(V0p,sub_vars,sub_val),x([1 3]),[Q(:)';QD(:)']),size(Q,1),[]);
h_X_val = reshape(msubs(subs(h_X,sub_vars,sub_val),x([1;3]),[Q(:)';QD(:)']),size(Q,1),[]);
h_Xp_val = reshape(msubs(subs(subs(h_X,x,xp),sub_vars,sub_val),x([1;3]),[Q(:)';QD(:)']),size(Q,1),[]);

Wval = reshape(msubs(subs(Wsol,sub_vars,sub_val),x([1;3]),[Q(:)';QD(:)']),size(Q,1),[]);
figure(1)
hold off
[cl,h]=contour(Q,QD,Wval,[1 1],'b');
clabel(cl,h);
title('W')
hold on
[cl,h]=contour(Q,QD,h_goal_val,[0 0],'k');
[cl,h]=contour(Q,QD,h_X_val,[0 0],'r');
[cl,h]=contour(Q,QD,h_Xp_val,[0 0],'r');
hold off

Vval = reshape(msubs(subs(Vsol,sub_vars,sub_val),x([1;3]),[Q(:)';QD(:)']),size(Q,1),[]);
Vval2 = reshape(msubs(subs(subs(Vsol,x,xp),sub_vars,sub_val),x([1;3]),[Q(:)';QD(:)']),size(Q,1),[]);

figure(ind*10+12)
hold off
[cl,h]=contour(Q,QD,Vval,[0 0],'b');
clabel(cl,h);
title('V')
hold on
[cl,h]=contour(Q,QD,Vval2,[0 0],'g');
[cl,h]=contour(Q,QD,h_goal_val,[0 0],'k');
[cl,h]=contour(Q,QD,h_X_val,[0 0],'r');
[cl,h]=contour(Q,QD,h_Xp_val,[0 0],'r');
hold off

%% sim test

sample_plot = true;

Vdotsol = sol.eval(Vdot);
dVdusol = clean(diff(Vdotsol,u));
% x0 = [-.55;0;.7;0];
% x0 = [0;0;.3;0];
% x0 = [-1.5;0;.9;0];
% x0 = [-.15;0;1;0];

if sample_plot
  [Q,QD] = meshgrid(linspace(-2,2,20),linspace(-2,2,20));
  Q = Q(:);
  QD = QD(:);
  VPVAL = Q*0;
else
  Q = .15;
  QD = 1;
end

for j=1:length(Q),
x0 = [Q(j);0;QD(j);0];

if double(subs(Vsol,[x;t],[x0;0])) > 0 && double(subs(h_X,x,x0)) > 0 && double(subs(subs(h_X,x,xp),x,x0)) > 0

% % x0 = [-.16;0;.6;0];
dt = 1e-2;
t_vec = 0:dt:T_unscaled;
x_vec = zeros(4,length(t_vec));
x_vec(:,1) = x0;
u_vec = zeros(1,length(t_vec));
% [vars,pows,coeffs] = decomp(f);

[vars,pows,coeffs] = decomp(dVdusol);
% [vars,pows,coeffs] = decomp(diff(V0p,u));

if time_varying
  inds = [1;4;3;5;2];
  assert(double(sum(vars(inds) - [x;t])) == 0)
else
  inds = [1;3;2;4];
  assert(double(sum(vars(inds) - x)) == 0)
end
pows = pows(:,inds);

for i=2:length(t_vec),
%   dVdu = double(subs(dVdusol,x,x_vec(:,i-1)));
if time_varying
  dVdu = coeffs*prod(repmat([x_vec(:,i-1)' t_vec(i)/T_unscaled],length(coeffs),1).^pows,2);
else
  dVdu = coeffs*prod(repmat([x_vec(:,i-1)'],length(coeffs),1).^pows,2);
end
  if dVdu > 0
    u_vec(i) = 1;
  else
    u_vec(i) = -1;
  end
  fi = [x_vec([3;4;1],i-1);u_vec(i)];
  x_vec(:,i) = x_vec(:,i-1) + dt*fi;
%   x_vec(:,i);
%   i
end
Vtraj = msubs(Vsol,[x;t],[x_vec;t_vec/T_unscaled]);
V0ptraj = msubs(V0p,[x;t],[x_vec;t_vec/T_unscaled]);
% figure(3)
% plot(t_vec,Vtraj,t_vec,V0ptraj);
% ylim([V0ptraj(1) 5])
% grid on



if sample_plot
  VPVAL(j) = max(V0ptraj);
  display(sprintf('%d:%f',j,VPVAL(j)))
else
  max(V0ptraj)
end

else
  VPVAL(j) = -inf;
end

end

if sample_plot
  figure(ind*10+12)
  hold on
  plot(Q(VPVAL>0),QD(VPVAL>0),'go')
  plot(Q(VPVAL<0),QD(VPVAL<0),'rx')
end
%% extract for saving
[dVduvars,dVdupows,dVducoeffs] = decomp(dVdusol);
if time_varying
  inds = [1;4;3;5;2];
  assert(double(sum(dVduvars(inds) - [x;t])) == 0)
else
  inds = [1;3;2;4];
  assert(double(sum(dVduvars(inds) - x)) == 0)
end
dVdupows = dVdupows(:,inds);


[Vpvars,Vppows,Vpcoeffs] = decomp(V0p);
inds = [1;3;2];
assert(double(sum(Vpvars(inds) - [q;v(1)])) == 0)
Vppows = [Vppows(:,inds) zeros(size(Vppows,1),1)];

[dVpdxvars,dVpdxpows,dVpdxcoeffs] = decomp(diff(V0p,x));
inds = [1;3;2];
assert(double(sum(dVpdxvars(inds) - [q;v(1)])) == 0)
dVpdxpows = [dVpdxpows(:,inds) zeros(size(dVpdxpows,1),1)];

save controller_V2 dVducoeffs dVdupows Vpcoeffs Vppows dVpdxcoeffs dVpdxpows

%%
save(sprintf('V%d',ind),'Vsol')