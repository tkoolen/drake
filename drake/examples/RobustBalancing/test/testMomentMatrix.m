function testMomentMatrix()
addpath_spotless;
h = addpathTemporary('..');

num_vars = 2;
for degree = 2 % TODO: 0, 1 doesn't work
  x = msspoly('x', num_vars);
  basis = grlex(monomials(x, 0 : 2 * degree));
  y = rand(length(basis), 1);
  
  [M, G] = momentMatrix(y, basis);
  moment_matrix_size = nchoosek(num_vars + degree, num_vars);
  sizecheck(M, [moment_matrix_size, moment_matrix_size]);

  p_coeffs = rand(moment_matrix_size, 1);
  p = p_coeffs' * grlex(monomials(x, 0 : degree));
  l_y_psquared = rieszFunctional(y, p * p, basis);
  l_y_psquared_back = p_coeffs' * M * p_coeffs;
  valuecheck(l_y_psquared, l_y_psquared_back);
end

end

function l_y = rieszFunctional(y, p, basis)
% From D. Henrion. Optimization on linear matrix inequalities for
% polynomial systems control. 
% http://homepages.laas.fr/henrion/Papers/henrion-grenoble14.pdf
%
% Given a sequence y = (y_\alpha)_{\alpha \in \mathbb{N}^n}, we define the
% Riesz linear functional l_y : \mathbb{R} \left[ x \right] \rightarrow
% \mathbb{R} such that l_y(x^\alpha) = y_\alpha for all \alpha \in
% \mathbb{N}^n.
%
% the ordering of y corresponds to basis

h = addpathTemporary('..');

if deg(p) == 0
  p_monomials = msspoly(1);
  coefficients = double(p);
else
  [vars, powers, coefficients] = decomp(p);
  p_monomials = recomp(vars, powers, eye(size(powers, 1)));
end
I = match_monomials(p_monomials, basis);

l_y = coefficients * y(I);

end
