function [M,G] = momentMatrix(y,b)
% 
%   y -- d-by-1  blah
%   b -- d-by-1  msspoly monomials.
    [x,p,C] = decomp(b);
    
    if any(sum(C,2) > 1), error('b must be monomials.'); end

    pow2 = unique(floor(p/2),'rows');
    phi = recomp(x,pow2,eye(size(pow2,1)));
    phi = grlex(phi);
    G = phi*phi';
    g = mss_s2v(G);
    
    I = match_monomials(g,b);
    if any(I == 0), error('Missing monomials.'); end
    m = y(I);
    M = mss_v2s(m);
    
end