function kordaPlot2d(model, sol, options)
clf;
legend_strings = {};
R_diag = options.R_diag;
xlim = [-R_diag(1) R_diag(1)];
ylim = [-R_diag(2) R_diag(2)];
hold on;

handles = [];

handles = [handles; contourSpotless(sol.g_X, sol.x(1), sol.x(2), xlim, ylim, [], [], 0, {'k'})];
legend_strings{end + 1} = 'g_X';

if isfield(sol, 'g_X_target')
  handles = [handles; contourSpotless(sol.g_X_target, sol.x(1), sol.x(2), xlim, ylim, [], [], 0, {'g'})];
  legend_strings{end + 1} = 'g_X^T';
end

handles = [handles; contourSpotless(sol.v_outer, sol.x(1), sol.x(2), xlim, ylim, [], [], 0, {'b'})];
legend_strings{end + 1} = 'v (outer)';

nbeta = length(options.betas);
colors = mat2cell(autumn(nbeta), ones(1, nbeta), 3);
inner_contour_handles = contourSpotless(-sol.vs_inner, sol.x(1), sol.x(2), xlim, ylim, [], [], zeros(nbeta, 1), colors);
% doesn't work because of Matlab bug:
%       handles = [handles; inner_contour_handles];
%       for i = 1 : nbeta
%         legend_strings{end + 1} = ['v (inner), \beta = ' num2str(options.betas(i))];
%       end

xlabel('x_1')
ylabel('x_2')
legend(gca(), handles, legend_strings);

while true
  x0 = zeros(model.num_states, 1);
  [x0(1), x0(2)] = ginput(1);
  lastwarn('');
  [~, x_traj] = ode45(@(t, x) msubs(sol.fbar, sol.x, x), [0 10], x0);
  if strcmp(lastwarn, '');
    plot(x_traj(:, 1), x_traj(:, 2), 'k-');
  end
end

hold off;
end
