function m = grlex(p)
    [x,p,M] = decomp(p);
    if any(sum(M==1,2) ~= 1)
        error('p must be monomials');
    end
    [~,I1] = sortrows(fliplr(p));
    [~,I2] = sort(sum(p(I1,:),2));
    m = recomp(x,p(I1(I2),:),speye(size(p,1)));
end