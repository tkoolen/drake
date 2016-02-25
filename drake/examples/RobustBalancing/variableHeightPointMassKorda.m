function variableHeightPointMassKorda()

g = 1;
step_max = .7;
step_time = 0.3;
z_nom = 1;
u_max = 1.1;
u_min = 0.1;

model = VariableHeightPointMass2D(g, z_nom, step_max, step_time, u_max, u_min);

options.R_diag = [2, .5 2, 2];
target_radius = 0.1;
%       target = @(x) 1 - (x / target_radius)' * (x / target_radius);
target = @(x) target_radius^2 - x' * x;

options.v_degree = 8;
options.w_v_degree = 10;
options.betas = 0.6; %[10, 0.6, 0.01, 0.001]; %0.1;
options.beta_outer_ind = 1;
sol = korda2014RegionOfAttractionAndController(model, target, options);

x = sol.x;

figure(1);
clf;
contourSpotless3D(subs(sol.v_outer, x(4), 0), [x(1); x(3); x(2)], 0, [options.R_diag(1); options.R_diag(3); options.R_diag(2)]);
title('outer approx');

figure(2);
clf;
contourSpotless3D(subs(-sol.vs_inner(end), x(4), 0), [x(1); x(3); x(2)], 0, [options.R_diag(1); options.R_diag(3); options.R_diag(2)]);
title('inner approx');

end
