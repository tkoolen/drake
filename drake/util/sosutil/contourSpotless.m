function h=contourSpotless(poly,x_var,y_var,x_range,y_range,const_vars,const_vals,contour_value,color,x_count,y_count,offset)
% function h = contourSpotless(poly,x_var,y_var,x_range,y_range,...
%               [const_vars],[const_vals],[contour_value],[color],[x_count],[y_count])
% Utility function for contour plotting of a spotless polynomial
% Optional arguments, marked [arg] above, can be ommited or set to []
%
% @input poly -- polynomials (nx1) to plot
% @input x_var -- variable (msspoly) for the x-axis
% @input y_var -- variable for the y-axis
% @input x_range -- [x_min x_max] range to plot
% @input y_range -- [y_min y_max] range to plot
% @input const_vars variables to hold constant in the plot
% @input const_vals corresponding constant values
% @input contour_value Set to plot only a specific value
% @input color Optional color argument
% @input x_count Number of x-values to sample, default 500
% @input y_count

poly = poly(:);
n_contour = length(poly);
if nargin < 6
  const_vars = [];
end
if nargin < 7
  const_vals = [];
end
if nargin < 8
  contour_value = [];
end
if nargin < 9
  color = [];
end
if nargin < 10
  x_count = 500;
end
if nargin < 11
  y_count = 500;
end

if nargin < 12
  offset = 0;
end


[X_VALS,Y_VALS] = meshgrid(linspace(x_range(1),x_range(2),x_count),linspace(y_range(1),y_range(2),y_count));
if ~isempty(const_vars)
  poly = subs(poly,const_vars,const_vals);
end

hold_value = ishold;

h = zeros(n_contour,1);

for i=1:n_contour
  POLY_VAL = reshape(msubs(poly(i),[x_var;y_var],[X_VALS(:)';Y_VALS(:)']),size(X_VALS,1),[]) + offset;
  
  if ~isempty(color) && ~isempty(contour_value)
    [cl,h(i)]=contour(X_VALS,Y_VALS,POLY_VAL,[contour_value(i) contour_value(i)]+offset,'Color',color{i});
  elseif ~isempty(color)
    [cl,h(i)]=contour(X_VALS,Y_VALS,POLY_VAL,color{i});
  elseif ~isempty(contour_value)
    [cl,h(i)]=contour(X_VALS,Y_VALS,POLY_VAL,[contour_value(i) contour_value(i)]+offset);
  else
    [cl,h(i)]=contour(X_VALS,Y_VALS,POLY_VAL);
  end
  clabel(cl);
  
  if i == 1
    hold on
  end
end

if ~hold_value
  hold off
end
end