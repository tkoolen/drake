function lipm2dKorda()

g = 1;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = 1; % set to 0 to get point foot model with no continuous inputs
model = LIPM2D(g, z_nom, step_max, step_time, cop_max);
options.R_diag = [1 1];
target_radius = 0.1;
%       target = @(x) 1 - (x / target_radius)' * (x / target_radius);
target = @(x) target_radius^2 - x' * x;

options.v_degree = 12;
options.w_v_degree = 12;
options.betas = [10, 0.6, 0.1, 0.01, 0.001]; %0.1;
options.beta_outer_ind = 2;
sol = korda2014RegionOfAttractionAndController(model, target, options);

figure(1);
kordaPlot2d(model, sol, options);
end
