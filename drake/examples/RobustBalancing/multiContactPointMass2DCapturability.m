function multiContactPointMass2DCapturability()
g = 1;
m = 1;
z_nom = 1;
step_max = 1;
step_time = 1;

% max_forces = [2, 2];
% normal_angles = [-deg2rad(70), deg2rad(70)];
% normals = [sin(normal_angles); cos(normal_angles)];
% mus = [0.1, 0.1];
% contact_points = [...
%   -0.5, 0.5;
%    0.0, 0.0];

% reduction to single point contact:
max_forces = 1.5;
normals = [0; 1];
mus = 3;
contact_points = [0; 0];

R_diag = [1, 1, 1, 1];

model = MultiContactPointMass2D(g, m, z_nom, step_max, step_time, max_forces, normals, mus, contact_points);

T = 2;
options.degree = 4;
options.do_backoff = false;
options.backoff_ratio = 1.04;

% R_diag = 2 * ones(1, model.num_states);

% goal_radius = 0.01;
% target = @(x) goal_radius^2 - x'*x;
target = [];

n = 0;
nStepCapturabilitySOS(model, T, R_diag, target, n, options)

end
