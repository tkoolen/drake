function contourSpotless3D(poly, vars, contour_value, R_diag)
% plot 3D contour (iso surface) of a polynomial
%
% @param poly polynomial to plot iso surface of
% @param vars 3-vector containing msspoly variables on the axes in the order x, y, z
% @param contour_value
% @param R_diag 3-vector containing radii of ellipsoid to show and 
% intersect with poly level set

sizecheck(vars, 3);

intersect_with_ellipsoid = true; % TODO: make option

colormap summer;

% uniform grid
grid_size = 50;
[x1s_grid, x2s_grid, x3s_grid] = meshgrid(...
  linspace(-R_diag(1), R_diag(1), grid_size),...
  linspace(-R_diag(2), R_diag(2), grid_size),...
  linspace(-R_diag(3), R_diag(3), grid_size));

% ellipsoid
sizecheck(R_diag, 3);
A = diag(1./(R_diag.^2));
h_X = 1 - vars'*A*vars;
h_Xs = full(msubs(h_X, vars, [x1s_grid(:)';x2s_grid(:)'; x3s_grid(:)']));
h_Xs = reshape(h_Xs, grid_size, grid_size, grid_size);

% evaluate polynomial on grid (always needed to use isonormals)
poly_vals_grid = reshape(full(msubs(poly, vars, [x1s_grid(:)';x2s_grid(:)'; x3s_grid(:)'])), grid_size, grid_size, grid_size);

x1s = x1s_grid;
x2s = x2s_grid;
x3s = x3s_grid;

if intersect_with_ellipsoid
  % project grid points outside of ellipsoid onto boundary of ellipsoid
  mask = h_Xs < 0;
  
  % x' * A * x = 1 - h_X
  % (s * x)' * A * (s * x) = 1
  % s^2 * (1 - h_X) = 1
  % s = sqrt(1 / (1 - h_X))
  s = sqrt(1 / (1 - h_Xs));
  x1s(mask) = s(mask) .* x1s_grid(mask);
  x2s(mask) = s(mask) .* x2s_grid(mask);
  x3s(mask) = s(mask) .* x3s_grid(mask);
end

% evaluate polynomial on region of interest (grid or ellipsoid)
poly_vals = reshape(full(msubs(poly, vars, [x1s(:)';x2s(:)'; x3s(:)'])), grid_size, grid_size, grid_size);

% actual zero-level set patch
c = colormap;
patch_poly = patch(isosurface(x1s, x2s, x3s, poly_vals, contour_value), 'FaceColor', c(1, :), 'EdgeColor', 'none', 'AmbientStrength', 0.7);

% patch that closes the region
patch_poly_caps = patch(isocaps(x1s,x2s,x3s,poly_vals,contour_value), 'FaceColor', 'interp', 'EdgeColor', 'none', 'AmbientStrength', 0.7);

% faint ellipsoid
patch(isosurface(x1s, x2s, x3s, h_Xs, contour_value), 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.03);

% use poly_val values to get better normals
isonormals(x1s_grid,x2s_grid,x3s_grid,poly_vals_grid, patch_poly);

if intersect_with_ellipsoid
  % use h_Xs normals for caps
  isonormals(x1s_grid,x2s_grid,x3s_grid,h_Xs, patch_poly_caps);
  
  % draw intersection between ellipsoid boundary and poly level set
  patch_poly_caps_boundary = patch(isocaps(x1s,x2s,x3s,poly_vals,contour_value), 'FaceColor', 'none');
  poly_vals_caps = full(msubs(poly, vars, patch_poly_caps_boundary.Vertices'));
  patch_poly_caps_boundary.FaceVertexCData = nan(length(poly_vals_caps), 3);
  tol = min(poly_vals_caps) + 1e-2 * max(poly_vals_caps) - min(poly_vals_caps);
  patch_poly_caps_boundary.FaceVertexCData(poly_vals_caps < tol, :) = 0;
%   patch_poly_caps_boundary.FaceVertexCData(abs(poly_vals_caps - contour_value) < tol, :) = 0;
  patch_poly_caps_boundary.EdgeColor = 'interp';
  patch_poly_caps_boundary.LineWidth = 2;
end

% formatting
daspect([1 1 1])
view(-37.5, 35)
light('Style', 'Local', 'Position', [2 0 -2]);
camlight
lighting gouraud
grid on;
box on;
ax = gca();
set(ax, 'BoxStyle', 'full');
ax.XTick = linspace(-R_diag(1), R_diag(1), 5);
ax.YTick = linspace(-R_diag(2), R_diag(2), 5);
ax.ZTick = linspace(-R_diag(3), R_diag(3), 5);
axis vis3d
% zoom(1.4)
% colorbar
end
