function [rho,sol] = lyapunovLevelSet(x,f,V,rho)
if nargin > 3
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  Vdot = diff(V,x)*f;
  
  [prog, Vdot_sos,mult,coeff] = spotless_add_sprocedure(prog, -Vdot, rho-V,x,4);
  prog = prog.withSOS(Vdot_sos);
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(coeff(1),solver,spot_options);
else
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  [prog,rho] = prog.newFree(1);
  Vdot = diff(V,x)*f;
  
  [prog, Vdot_sos] = spotless_add_sprocedure(prog, -(1+x'*x)^2*(V-rho) - Vdot, Vdot,x,0);
  prog = prog.withSOS(Vdot_sos);
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_sedumi;
  sol = prog.minimize(-rho,solver,spot_options);
  rho = sol.eval(rho);
end
end