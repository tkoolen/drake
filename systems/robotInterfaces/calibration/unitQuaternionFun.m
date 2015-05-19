function [g,dg] = unitQuaternionFun(quat)
  quat = reshape(quat,4,[]);
  num_quat = size(quat,2);
  g = sum(quat.^2,1)';
  iGfun = reshape(bsxfun(@times,ones(4,1),1:num_quat),[],1);
  jGvar = (1:4*num_quat)';
  Gval = 2*reshape(quat,[],1);
  dg = sparse(iGfun,jGvar,Gval,num_quat,4*num_quat);
end