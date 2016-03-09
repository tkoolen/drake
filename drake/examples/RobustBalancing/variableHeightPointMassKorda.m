function variableHeightPointMassKorda()

g = 1;
step_max = .7;
step_time = 0.3;
z_nom = 1;
u_max = 1.1;
u_min = 0.1;

model = VariableHeightPointMass2D(g, z_nom, step_max, step_time, u_max, u_min);

% options.R_diag = [2, .5 2, 2];
% % target_radius = 0.1;
% % target = @(x) 1 - (x / target_radius)' * (x / target_radius);
% % target = @(x) target_radius^2 - x' * x;
% target = zeros(model.num_states, 1);
% 
% options.v_degree = 8;
% options.w_v_degree = 8;
% options.betas = 0.0001; %[10, 0.6, 0.01, 0.001]; %0.1;
% options.beta_outer_ind = 1;
% sol = korda2014RegionOfAttractionAndController(model, target, options);
% 
% x = sol.x;
% 
% figure(1);
% clf;
% contourSpotless3D(subs(sol.v_outer, x(4), 0), [x(1); x(3); x(2)], 0, [options.R_diag(1); options.R_diag(3); options.R_diag(2)]);
% colorbar;
% title('outer approx');
% 
% figure(2);
% clf;
% contourSpotless3D(subs(-sol.vs_inner(end), x(4), 0), [x(1); x(3); x(2)], 0, [options.R_diag(1); options.R_diag(3); options.R_diag(2)]);
% colorbar;
% title('inner approx');
% 
% x0 = [0.1, 0.1, 0.1, 0.1];
% [~, x_traj] = ode45(@(t, x) msubs(sol.fbar, sol.x, x), [0 10], x0);

options.split_inputs = false;
options.beta = 1;
options.R_diag = [0.5, .5, 1, 1];
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
xlabel('x'); ylabel('z');
hold on;

while true
  x0 = zeros(model.num_states, 1);
  [x0(1), x0(2)] = ginput(1);
  x0(3) = -1.1 * x0(1);
  lastwarn('');
  T = 5;
  [~, x_traj] = ode4(@(t, x) sol.fbar(x), 0 : 0.01 : T, x0);
  if strcmp(lastwarn, '');
    plot(x_traj(:, 1), x_traj(:, 2), 'k-');
    axis([xlim, ylim]);
  end
end

end
