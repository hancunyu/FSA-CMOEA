function [Mate1,Mate2,Mate3] = opfpim3e_neighbor_pairing(MatingPop,Zmin)
%OPFPIU_NEIGHBOR_PAIRING Scale-invariant objective-direction pairing.
%
% Objective differences from the running ideal point are scaled per
% objective before cosine similarity is computed. This preserves the
% original neighbor-pairing logic while removing engineering-unit bias.

    Objs = MatingPop.objs;
    [Num,M] = size(Objs);
    Objs = Objs-repmat(Zmin,Num,1);

    % Scale each objective relative to the current ideal point.
    scale = max(Objs,[],1);
    threshold = 64*eps.*max(1,abs(scale));
    scale(scale<=threshold) = 1;
    Objs = Objs./repmat(scale,Num,1);
    Objs(~isfinite(Objs)) = 0;

    norms = sqrt(sum(Objs.^2,2));
    zeroMask = norms<=eps;
    norms(zeroMask) = 1;
    Objs = Objs./repmat(norms,1,M);
    if any(zeroMask)
        Objs(zeroMask,:) = 0;
    end

    CosV = Objs*Objs'-3*eye(Num);
    [~,SInd] = sort(-CosV,2);
    Nr = min(10,max(2,Num-1));
    Neighbor = SInd(:,1:Nr);

    Mate1 = MatingPop;
    P = ones(Num,2);
    for i = 1:Num
        if Nr>=2
            P(i,:) = Neighbor(i,randsample(Nr,2));
        else
            P(i,:) = Neighbor(i,[1,1]);
        end
    end
    Mate2 = MatingPop(P(:,1));
    Mate3 = MatingPop(P(:,2));
end





