clear all
N = 10;
T = 5;
duration = [T T];

% Trajectory optimization
g = 10;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = .1; % set to 0 to get point foot model with no continuous inputs
fz_range = .3;
fx_range = .3;
inertia_ratio = .6^2/2;
R_diag = [1 cop_max*3*sqrt(g/z_nom) 1 1 2 2];

model = TransformedFull2DModel(g, inertia_ratio, z_nom, step_max, step_time, fz_range, fx_range, cop_max);
A = diag(1./(R_diag.^2));

plant = NStepCapturabilityPlant(model,false);
plant = plant.setInputLimits(-ones(3,1),ones(3,1));

[X1,X2] = meshgrid(linspace(-1,1,30),linspace(-1,1,30));
X1 = X1(:);
X2 = X2(:);
info_vec = zeros(size(X1));
cost_vec = zeros(size(X1));
for i=1:length(X1),
%   x0 = [0;.55;0;0;0;0];
  x0 = [X1(i);X2(i);0;0;0;0];
  
  if x0'*A*x0 < 1
  
  traj_opt = DircolTrajectoryOptimization(plant,N,duration);
  
  
  traj_opt = traj_opt.setSolverOptions('snopt','print','snopt.out');
  traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',50);
  
  % traj_opt = traj_opt.addRunningCost(
  
  init_constraint = ConstantConstraint(x0);
  traj_opt = traj_opt.addStateConstraint(init_constraint,1);
  
  final_cost = QuadraticConstraint(-inf,inf,10*eye(6),zeros(6,1));
  traj_opt = traj_opt.addCost(final_cost,traj_opt.x_inds(:,end));
  
  state_ball = QuadraticConstraint(0,.5,A,zeros(6,1));
  traj_opt = traj_opt.addStateConstraint(state_ball,2:N);
  
  t_init = linspace(0,T,N);
  traj_init.x = PPTrajectory(foh(t_init,repmat(x0,1,N)));
  
  tic
  [xtraj,utraj,z,F,info] = traj_opt.solveTraj(t_init,traj_init);
  toc
  [info F(1)]
  
  info_vec(i) = info;
  cost_vec(i) = F(1);
  else
    info_vec(i) = inf;
    cost_vec(i) = inf;
  end
  i
end
%%
v = NStepCapturabilityVisualizer(plant);
v = v.setInputFrame(plant.getStateFrame);
v.playback_speed = T;
figure(25)
v.playback(xtraj);

I = info_vec < 10 & cost_vec < 1;