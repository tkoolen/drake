function [prog, eqn, mult, coefmult] = spotless_add_sprocedure(prog, eqn, h, vars, degree,sos_option)
%SPOTLESS_ADD_SPROCEDURE Summary of this function goes here

% eqn_deg = full(deg(eqn,vars));
% if ~even(eqn_deg)
%   eqn_deg = eqn_deg + 1;
% end
% degree = max(2,min(degree, 2*floor((eqn_deg - deg(h))/2)));

if nargin < 6
  sos_option = 1;
end

original_deg = even_degree(eqn,vars);

if isempty(degree)
  degree = zeros(length(h),1);
  for i = 1:length(h)
    degree(i) = original_deg - even_degree(h(i),vars);
  end
else
  if length(degree) == 1
    degree = repmat(degree,length(h),1);
  end
  
end

mult = msspoly;
coefmult = msspoly;
for i = 1 : length(h)
  switch sos_option
    case 1
      [prog,mult_i,coefmult_i] = prog.newSOSPoly(monomials(vars,0:degree));
    case 2
      [prog,mult_i,coefmult_i] = prog.newSDSOSPoly(monomials(vars,0:degree));
    case 3
      [prog,mult_i,coefmult_i] = prog.newDSOSPoly(monomials(vars,0:degree));
  end  
  
  eqn = eqn - h(i) * mult_i;
  mult = [mult; mult_i]; %#ok<AGROW>
  coefmult = [coefmult; coefmult_i]; %#ok<AGROW>
  
  display(sprintf('S-proc ineq. SOS deg: %d, h deg: %d, mult deg: %d',original_deg, deg(h(i)), degree(i)))
  if original_deg ~= even_degree(h(i),vars) + degree(i);
    warning('S-procedure degree mismatch')
  end
end
end

