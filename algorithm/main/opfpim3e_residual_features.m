function F = opfpim3e_residual_features(PopCon,profile)
%OPFPIU_RESIDUAL_FEATURES Physically grouped constraint-violation features.
%
% The first two columns are always CVsum and CVmax. Additional columns are
% group-wise sums and maxima when a benchmark constraint map is available.

    V = max(0,PopCon);
    if isempty(V)
        F = zeros(0,2);
        return;
    end

    Fsum = sum(V,2);
    Fmax = max(V,[],2);
    F = [Fsum,Fmax];

    groups = profile.constraintGroups;
    for k = 1:numel(groups)
        idx = groups{k};
        idx = idx(idx>=1 & idx<=size(V,2));
        if isempty(idx)
            F = [F,zeros(size(V,1),2)]; %#ok<AGROW>
        else
            F = [F,sum(V(:,idx),2),max(V(:,idx),[],2)]; %#ok<AGROW>
        end
    end

    if isempty(groups) && size(V,2)>0
        half = floor(size(V,2)/2);
        if half > 0
            F = [F,sum(V(:,1:half),2),...
                sum(V(:,half+1:end),2)]; %#ok<AGROW>
        end
    end
end





