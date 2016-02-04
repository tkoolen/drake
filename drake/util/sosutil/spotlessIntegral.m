function cost = spotlessIntegral(prog,poly,sphere_vars,A_diag,box_vars,box_lims)
% l = spotlessIntegral(poly,sphere_vars,A,box_vars,box_lims)

[vars,alphas,coeff] = decomp(poly, prog.freeVar);

% find sphere_vars
n_sphere = length(sphere_vars);

if n_sphere > 0
  valuecheck(length(A_diag),n_sphere);
end

sphere_inds = [];
for i=1:length(sphere_vars)
  for j=1:length(vars)
    if isequal(vars(j),sphere_vars(i))
      sphere_inds(end+1) = j;
      break;
    end
  end
end

% find box
n_box = length(box_vars);
if n_box > 0
  sizecheck(box_lims,[n_box 2]);
end

box_inds = [];
for i=1:length(box_vars)
  for j=1:length(vars)
    if isequal(vars(j),box_vars(i))
      box_inds(end+1) = j;
      break;
    end
  end
end

assert(length([sphere_inds;box_inds]) == length(unique([sphere_inds;box_inds])));
assert(length([sphere_inds;box_inds]) == length(vars))

sphere_alphas = alphas(:,sphere_inds);
box_alphas = alphas(:,box_inds);

% eliminate box variables first
% /int x^alpha = 1/(alpha+1)*(x_max^(alpha+1) - x_min^(alpha+1))
for i=1:length(box_inds)
  pows = box_alphas(:,box_inds(i)) + 1;
  vals = 1./pows .* (box_lims(i,2).^pows - box_lims(i,1).^pows);
  coeff = coeff.*vals';
end

if n_sphere > 0
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