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
  
  display(sprintf('S-proc ineq. SOS deg: %d, h deg: %d, mult deg: %d',full(deg(eqn,vars)), deg(h(i)), degree))
end
end

