load V0_TransformedFull2DModel
options.zero_origin = true;
u_bnd = boundController(model,uCoeff,uMomentBasis,T,R_diag,options);

%% Get an initial quadratic Lyapunov candidate
x = msspoly('x',model.num_states);
t = msspoly('t',1);
f = model.dynamics(t,x,2*u_bnd-1);

prog = spotsosprog();
[prog,R] = prog.newPSD(model.num_states);
[prog,rho] = prog.newPos(1);

% calculate hessian of Vdot
V = x'*(R + rho*eye(model.num_states))*x;
Vdot = diff(V,x)*f;
H=.5*subs(diff(diff(Vdot,x)',x),[t;x],zeros(model.num_states+1,1));
prog = prog.withEqs(R(1) - 1);

prog = prog.withPSD(-H - rho*eye(model.num_states));

spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = true;
solver = @spot_mosek;
sol = prog.minimize(-rho,solver,spot_options);

V0 = sol.eval(V);
%%
[ V_inner] = quadraticLyapunovAlternations(x,f,V0*100)

%%
u_bnd = 2*u_bnd - 1;
save('V0_TransformedFull2DModel_inner','Vsol','model','T','R_diag','u_bnd','V_inner')
%%
load V0_TransformedFull2DModel_inner

n = 1;

options.degree = 4;
options.scale = 1;
options.control_design = true;
options.korda_control_design = true;
options.beta = 1;
options.infinite_time = false;
options.free_final_time = false;
options.V0 = V_inner;

goal_radius = 0.05;
target = @(x) goal_radius^2 - x'*x;

[Vsol,Wsol] = nStepCapturabilitySOS(model, T, R_diag, target, n, options);

%%
% options = struct();
% x=msspoly('x', model.num_states); % state
% goal_radius = 0.05;
% target = @(x) goal_radius^2 - x'*x;
% options.target = target;
% options.degree = 4;
% options.Wlvl = .2;
% 
% [Vsol2,Vlvl] = altLyapunovInnerApproximation(model,Vsol,T,2*u_bnd-1,R_diag,options);
% 
% % options.beta = 1;
% 
% % [Vsol2,Wsol2] = innerApproximation(model,u_sol,R_diag,options.target,options);