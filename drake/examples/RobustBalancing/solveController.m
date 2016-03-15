function [u,coeff] = solveController(M1,M2,G,condTh)
if nargin < 4, condTh = 1e6; end
% Now analyze M1.
if condTh < Inf
  %       u = G(1,:)*pinv(M1,1/condTh)*M2(:,1);
  %       return;
  [U,S,V] = svd(M1);
  s = diag(S);
  keep = find((s(1)./s) < condTh,1,'last');
  
  coeff = U\M2(:,1);
  coeff = diag(1./s(1:keep))*coeff(1:keep,:);
  coeff = V(:,1:keep)*coeff;
  u = G(1,:)*coeff;
else
  coeff = (M1\M2(:,1));
  u = G(1,:)*coeff;
end
end
