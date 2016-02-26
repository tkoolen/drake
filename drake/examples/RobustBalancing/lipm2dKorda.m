function lipm2dKorda()

g = 1;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = 0.9; % set to 0 to get point foot model with no continuous inputs
model = LIPM2D(g, z_nom, step_max, step_time, cop_max);


% options.R_diag = [1 1];
% target_radius = 0.1;
% target = @(x) 1 - (x / target_radius)' * (x / target_radius);
% target = @(x) target_radius^2 - x' * x;
% target = zeros(model.num_states, 1);

% options.v_degree = 12;
% options.w_v_degree = 12;
% options.betas = [10, 0.6, 0.1, 0.01, 0.001]; %0.1;
% options.beta_outer_ind = 2;
% sol = korda2014RegionOfAttractionAndController(model, target, options);
% 
% figure(1);
% kordaPlot2d(model, sol, options);

options.split_inputs = true;
options.beta = 1;
options.R_diag = [1, 1];
options.ubar = 1;
options.l_x = @(x) x' * x;
options.l_u = @(x) 1e-1;
options.degree = 10;
options.M = 1.01 * (options.l_x([1; 0]) + options.l_u([1; 0]) * options.ubar) / options.beta;
sol = korda2015(model, options);

figure(1);
clf;
R_diag = options.R_diag;
xlim = [-R_diag(1) R_diag(1)];
ylim = [-R_diag(2) R_diag(2)];
axis([xlim, ylim]);
hold on;

while true
  x0 = zeros(model.num_states, 1);
  [x0(1), x0(2)] = ginput(1);
  lastwarn('');
  T = 10;
  [~, x_traj] = ode45(@(t, x) sol.fbar(x), [0 T], x0);
  if strcmp(lastwarn, '');
    plot(x_traj(:, 1), x_traj(:, 2), 'k-');
    axis([xlim, ylim]);
  end
end

end
