function Del = opfpim3e_truncation(PopObj,K)
%OPFPIU_TRUNCATION Scale-invariant SPEA2 truncation.

    n = size(PopObj,1);
    Del = false(1,n);
    if K<=0 || n==0
        return;
    end

    Geometry = opfpim3e_scale_geometry(PopObj);
    Distance = pdist2(Geometry,Geometry);
    Distance(logical(eye(n))) = inf;

    while sum(Del)<min(K,n)
        Remain = find(~Del);
        if isempty(Remain)
            break;
        end
        Temp = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end





