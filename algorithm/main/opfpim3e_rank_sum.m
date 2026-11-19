function score = opfpim3e_rank_sum(X)
% Parameter-free rank aggregation; ties receive stable ordinal ranks.
    n=size(X,1); score=zeros(n,1);
    for j=1:size(X,2)
        [~,ord]=sort(X(:,j),'ascend');
        r=zeros(n,1); r(ord)=1:n;
        score=score+r/max(1,n);
    end
end





