function [prog, eqn, mult, coefmult ] = spotless_add_eq_sprocedure(prog, eqn, g, vars, degree)
%SPOTLESS_ADD_EQ_SPROCEDURE Summary of this function goes here

% eqn_deg = full(deg(eqn,vars));
% if ~even(eqn_deg)
%   eqn_deg = eqn_deg + 1;
% end
% degree = min(degree, eqn_deg - deg(g));

[prog,mult,coefmult] = prog.newFreePoly(monomials(vars,0:degree));
eqn = eqn - g*mult;

display(sprintf('S-proc eq. SOS deg: %d, g deg: %d, mult deg: %d',full(deg(eqn,vars)), deg(g), degree))

end

