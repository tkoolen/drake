function [ V,u ] = quadraticControlLyapunovAlternations(x_mss,u_mss,f,V0)
N = 1;
nX = length(x_mss);
V = V0;
for i=1:1,
  prog = spotsosprog;
  prog = prog.withIndeterminate(x_mss);
  [prog,gamma] = prog.newPos(1);
  [prog,u] = prog.newFreePoly(monomials(x_mss,1:3),3);
  Vdot = diff(V,x_mss)*subs(f,u_mss,u);

  % not obeying input limits
  
  [prog, Vdot_sos,mult,coeff] = spotless_add_sprocedure(prog, -Vdot, 1-V,x_mss,4);
  prog = prog.withSOS(Vdot_sos - gamma*(1+x_mss'*x_mss)^3);
  for i=1:length(u),
    [prog, u_sos,upmult{i}] = spotless_add_sprocedure(prog, 1-u(i), 1-V,x_mss,2);
    prog = prog.withSOS(u_sos);
    [prog, u_sos,ummult{i}] = spotless_add_sprocedure(prog, u(i)+1, 1-V,x_mss,2);
    prog = prog.withSOS(u_sos);
  end
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(-gamma,solver,spot_options);
  u = sol.eval(u)
  
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
  
%   keyboard
  
  solver = @spot_mosek;
  prog = prog.withPos(1.03*sol.eval(trace(Q)) - trace(Q));
%   sol = prog.minimize(gamma,solver,spot_options);
  
  V = sol.eval(V);

  %   keyboard
  
end
end

