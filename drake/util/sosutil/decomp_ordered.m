function [pows,coeffs] = decomp_ordered(p,ordered_vars,v)
 % Run spotless decomp(p,v) but return arguments in order given by x
 if nargin > 2
  [vars,pows,coeffs] = decomp(p,v);
 else
   [vars,pows,coeffs] = decomp(p);
 end
 n_vars = length(vars);
 n_ordered = length(ordered_vars);
 
 inds = (n_vars+1)*ones(n_ordered,1); % initialize to be n_vars+1
 for i=1:n_ordered,
   for j=1:n_vars,
     if isequal(ordered_vars(i), vars(j))
       inds(i) = j;
       break;
     end
   end
 end
 
 % add an extra zero variable and pows for elements in ordered_vars that
 % are not in vars
 vars = [vars;0];
 pows = [pows zeros(size(pows,1),1)];
 
 vars = vars(inds);
 pows = pows(:,inds);
end