function degree = even_degree(p,vars)
  % deg = even_degree(p,vars)
  % calculate the degree of polynomial p, rounded up to the nearest even
  % number
  degree = full(2*ceil(deg(p,vars)/2));
end