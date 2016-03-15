function transformedLipm2DNstepCapturability(n)
if nargin < 1
  n = 0;
end
g = 10;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = .1; % set to 0 to get point foot model with no continuous inputs

model = TransformedLIPM2D(g, z_nom, step_max, step_time, cop_max);
R_diag = [1 cop_max*2*sqrt(g/z_nom)];
if n > 0
  T = step_time;
else
  T = 2;
end
options.degree = 6;
options.scale = 1;
options.control_design = true;
options.korda_control_design = true;
options.beta = 1;
options.infinite_time = true;
options.free_final_time = false;

% radius of ball around the origin used as goal for 0-step capturability
goal_radius = 0.01;
target = @(x) goal_radius^2 - x'*x;
target = [];

[Vsol,Wsol] = nStepCapturabilitySOS(model, T, R_diag, target, n, options);

% keyboard
%%

% 
% 
% options.beta = [10 1 .1 .01 .001];
% % options.beta = .1;
% V_inner = innerApproximation(model,u_sol,R_diag,target,options);
% 
% x = msspoly('x',model.num_states);
% data{1}.u_sol = u_sol;
% data{1}.T = T;
% plant = HybridCapturabilityPlant(model,data);
% 
% figure(1)
% kordaPlot2d(model,plant, x,1 - x'*diag(1./R_diag.^2)*x,target(x), Vsol, sum(V_inner), R_diag);
% 
% % keyboard
end
