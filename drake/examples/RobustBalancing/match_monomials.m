function I = match_monomials(g,b)
%
% I = match_monomials(g,b)
%
%  g -- n-by-m msspoly of monomials.
%  b -- k-by-1 msspoly of monomials.
%
%  I -- n-by-m msspoly where g(i) == b(I(i)).
%                      or I(i) = 0 if no match exists.
%
    
    [xb,pb,Cb] = decomp(b);
    [xg,pg,Cg] = decomp(g);
    
    if ~all(Cb(:) >= 0) || ~all(Cg(:) >= 0) || ...
            any(sum(Cb,2) > 1) || any(sum(Cg,2) > 1)
        error('Inputs must be vectors of monomials.');
    end
    
    % First, extend both to have the full set of variables.
    xall = decomp([xb;xg]);
    
    pball = extend_rep(xall,xb,pb);
    pgall = extend_rep(xall,xg,pg);
    

    %  Now, for each row in pgall, determine where it appears in pball.
    mall = sparse(ones(size(pb,1),size(pg,1)));
    for i = 1:length(xall)
        mall = mall.*(repmat(pball(:,i),1,size(pgall,1))==repmat(pgall(:,i)',size(pball,1),1));
    end
    
    [i,j] = ind2sub(size(mall),find(mall));
    % i refers to which monomials appearing in pb.
    
    brows = max(Cb.*repmat((1:size(Cb,1))',1,size(pb,1)));
    
    I = zeros(size(pg,1),1);
    I(j) = brows(i);
    I = Cg*I;
    I = reshape(I,size(g));

    function [pall] = extend_rep(xall,xa,pa)
        m = match(xall,xa);
        pall = sparse([],[],[],size(pa,1),length(xall));
        pall(:,m) = pa;
    end
end