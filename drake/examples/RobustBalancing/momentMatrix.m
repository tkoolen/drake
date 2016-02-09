function [M,G] = momentMatrix(y, basis)
% 
%   y -- d-by-1 moments ordered to correspond to grlex(basis)
%   basis -- d-by-1  msspoly monomial basis.
    [x,p,C] = decomp(basis);
    
    if any(sum(C,2) > 1), error('basis must be monomials.'); end

    pow2 = unique(floor(p/2),'rows');
    phi = recomp(x,pow2,eye(size(pow2,1)));
    phi = grlex(phi);
    G = phi*phi';
    g = mss_s2v(G);
    
    I = match_monomials(g,basis);
    if any(I == 0), error('Missing monomials.'); end
    m = y(I);
    M = mss_v2s(m);
    
end