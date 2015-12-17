function u = solveController(M1,M2,G,condTh)
    if nargin < 4, condTh = Inf; end
    % Now analyze M1.
    if condTh < Inf
        [U,S,V] = svd(M1);
        s = diag(S);
        keep = find((s(1)./s) < condTh,1,'last');

        v = U\M2(:,1);
        v = diag(1./s(1:keep))*v(1:keep,:);
        v = V(:,1:keep)*v;
        u = G(1,:)*v;
    else
        u = G(1,:)*(M1\M2(:,1));
    end
end