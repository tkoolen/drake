% Find V(t,x), rho(t), u(t,x) and s(x)
%  s.t.
%   V - rho => Vdot<= rhodot
%       sig1*(V - rho) - Vdot + rhodot >= 0
%   r(x,s(x)) = xp and V(T,x) <= rho(T) ==> Vp(0,xp) <= 1
%     p*(r(x,s(x)) - xp) + sig2*(V(T,x) - rho(T)) + 1 - Vp(0,xp) >= 0
%
%   rho(0) = 1;
%
%  Alternation A: sig1, sig2, u, s
%
%  Alternation B: V, rho, p

% Find V(x) st
%  V(x) <= 1 ==> Vp(r(x,s)) <= 1

min_s [x+As]'*S*[x+As]= x'Sx + s'A'SAs + 2x'SAs
s = (A'SA)^-1 * A'Sx
xp = x + A(A'SA)^-1 * A'Sx
xp = Tx
Vp(xp) = x'T'Q*T*x

function [ V,u ] = quadraticControlAlternationsWithReset(x_mss,u_mss,s_mss,f,r,V0,T,Vp)
u_degree = 3;
nX = length(x_mss);
nS = length(s_mss);
V = V0;
prog = spotsosprog;
prog = prog.withIndeterminate(x_mss);
[prog,t] = prog.newIndeterminate('t',1);
[prog,xp] = prog.newIndeterminate('xp',nX);
[prog,gamma] = prog.newPos(1);

[prog,u] = prog.newFreePoly(monomials([t;x_mss],1:u_degree),length(u_mss));
% [prog,s] = prog.newFreePoly(monomials(x_mss,1:u_degree),nS);
s = -x_mss(2);
Vdot = diff(V,x_mss)*subs(f,u_mss,u) + diff(V,t);

% satisfy input limits
[prog, Vdot_sos,mult1] = spotless_add_sprocedure(prog, -Vdot, 1-V,[t;x_mss],4);
prog = prog.withSOS(Vdot_sos - gamma*(1+x_mss'*x_mss)^3);
for i=1:length(u),
  [prog, u_sos,upmult{i}] = spotless_add_sprocedure(prog, 1-u(i), 1-V,[t;x_mss],2);
  prog = prog.withSOS(u_sos);
  [prog, u_sos,ummult{i}] = spotless_add_sprocedure(prog, u(i)+1, 1-V,[t;x_mss],2);
  prog = prog.withSOS(u_sos);
end

[prog, VT_sos,mult2] = spotless_add_sprocedure(prog, -Vdot, 1-subs(V,t,T),[x_mss;xp],4);
[prog, VT_sos,eqmult] = spotless_add_eq_sprocedure(prog, VT_sos, subs(r,s_mss,s) - xp,[x_mss;xp],4);
prog = prog.withSOS(VT_sos);
% for i=1:length(s),
%   [prog, s_sos,spmult{i}] = spotless_add_sprocedure(prog, 1-s(i), 1-subs(V,t,T),[t;x_mss],2);
%   prog = prog.withSOS(s_sos);
%   [prog, s_sos,smmult{i}] = spotless_add_sprocedure(prog, s(i)+1, 1-subs(V,t,T),[t;x_mss],2);
%   prog = prog.withSOS(s_sos);
% end

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(-gamma,solver,spot_options);
u = sol.eval(u)

keyboard

for i=1:length(u),
  upmult{i} = sol.eval(upmult{i});
  ummult{i} = sol.eval(ummult{i});
end
%   keyboard

mult = sol.eval(mult);

prog = spotsosprog;
prog = prog.withIndeterminate(x_mss);
[prog,gamma] = prog.newPos(1);
[prog,Q] = prog.newPSD(nX);
V = x_mss'*Q*x_mss;
Vdot = diff(V,x_mss)*subs(f,u_mss,u);

prog = prog.withSOS(-Vdot + mult*(V-1) + (1+x_mss'*x_mss)^3*1e-6);

for i=1:length(u),
  prog = prog.withSOS(upmult{i}*(V-1) + 1 - u(i));
  prog = prog.withSOS(ummult{i}*(V-1) + 1 + u(i));
end

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(trace(Q),solver,spot_options);

V = sol.eval(V);
end