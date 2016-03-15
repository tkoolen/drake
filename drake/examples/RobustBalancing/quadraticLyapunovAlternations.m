function [ V ] = quadraticLyapunovAlternations(x,f,V0)
N = 1;
nX = length(x);
V = V0;
for i=1:5,
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  [prog,gamma] = prog.newPos(1);
  Vdot = diff(V,x)*f;
  
  [prog, Vdot_sos,mult,coeff] = spotless_add_sprocedure(prog, -Vdot, 1-V,x,4);
  prog = prog.withSOS(Vdot_sos - gamma*(1+x'*x)^3);
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = true;
  spot_options.do_fr = true;
  solver = @spot_mosek;
  sol = prog.minimize(-gamma,solver,spot_options);
  
%   keyboard
  
  mult = sol.eval(mult);
  
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  [prog,gamma] = prog.newPos(1);
  [prog,Q] = prog.newPSD(nX);
  V = x'*Q*x;
  Vdot = diff(V,x)*f;
  
  prog = prog.withSOS(-Vdot + mult*(V-1) + (1+x'*x)^3*1e-6);
  
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

