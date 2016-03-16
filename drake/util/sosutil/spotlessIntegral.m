function cost = spotlessIntegral(prog,poly,sphere_vars,A_diag,box_vars,box_lims)
% l = spotlessIntegral(poly,sphere_vars,A,box_vars,box_lims)

[vars,alphas,coeff] = decomp(poly, prog.freeVar);

% find sphere_vars
n_sphere = length(sphere_vars);

if n_sphere > 0
  valuecheck(length(A_diag),n_sphere);
end

sphere_inds = [];
sphere_var_inds = [];
sphere_var_inds_c = [];
for i=1:length(sphere_vars)
  is_var = false;
  for j=1:length(vars)
    if isequal(vars(j),sphere_vars(i))
      sphere_inds(end+1) = j;
      sphere_var_inds(end+1) = i;
      is_var = true;
      break;
    end
  end
  if ~is_var
    sphere_var_inds_c(end+1) = i;
  end
end

% find box
n_box = length(box_vars);
if n_box > 0
  sizecheck(box_lims,[n_box 2]);
end

box_inds = [];
box_var_inds = [];
box_var_inds_c = [];
for i=1:length(box_vars)
  is_var = false;
  for j=1:length(vars)
    if isequal(vars(j),box_vars(i))
      box_inds(end+1) = j;
      box_var_inds(end+1) = i;
      is_var = true;
      break;
    end
  end
  if ~is_var
    box_var_inds_c(end+1) = i;
  end
end

assert(length([sphere_inds box_inds]) == length(unique([sphere_inds box_inds])));
assert(length([sphere_inds box_inds]) == length(vars))

sphere_alphas = zeros(size(alphas,1),length(sphere_vars));
sphere_alphas(:,sphere_var_inds) = alphas(:,sphere_inds);

box_alphas = zeros(size(alphas,1),length(box_vars));
box_alphas(:,box_var_inds) = alphas(:,box_inds);

% eliminate box variables first
% /int x^alpha = 1/(alpha+1)*(x_max^(alpha+1) - x_min^(alpha+1))
for i=1:length(box_vars)
  pows = box_alphas(:,i) + 1;
  vals = 1./pows .* (box_lims(i,2).^pows - box_lims(i,1).^pows);
  coeff = coeff.*vals';
end

if n_sphere > 0
  if isempty(alphas)
    sphere_alphas = zeros(size(A_diag));
  end
  
  sphere_betas = 0.5*(sphere_alphas + 1);
  Ra = (1.^(sum(sphere_alphas,2) + n_sphere))./(sum(sphere_alphas,2) + n_sphere);
  IS = 2*prod(gamma(sphere_betas),2)./(gamma(sum(sphere_betas,2)));
  l = Ra.*IS;
  alphaszero = (mod(sphere_alphas,2) ~= 0);
  alphaszero = any(alphaszero,2);
  l(alphaszero) = 0;
  
  l = l.*prod(repmat(A_diag,size(sphere_alphas,1),1).^(sphere_alphas+1),2);
  
  cost = coeff*l;
else
  cost = sum(coeff);
end
end