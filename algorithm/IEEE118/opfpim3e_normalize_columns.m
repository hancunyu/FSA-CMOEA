function Y = opfpim3e_normalize_columns(X)
    if isempty(X), Y=X; return; end
    lo=min(X,[],1); hi=max(X,[],1); scale=hi-lo;
    scale(scale<=eps)=1;
    Y=(X-repmat(lo,size(X,1),1))./repmat(scale,size(X,1),1);
end





