function Fitness = opfpim3e_calfitness(PopObj,PopCon)
%OPFPIU_CALFITNESS SPEA2-style fitness with scale-invariant density.
%
% Dominance is computed from the original objective values and constraint
% violations. Only the Euclidean geometry used by the density term is
% min-max scaled, so engineering units cannot dominate diversity pressure.

    N = size(PopObj,1);
    if nargin < 2 || isempty(PopCon)
        CV = zeros(N,1);
    else
        CV = sum(max(0,PopCon),2);
    end

    Dominate = false(N);
    for i = 1:N-1
        for j = i+1:N
            if CV(i) < CV(j)
                Dominate(i,j) = true;
            elseif CV(i) > CV(j)
                Dominate(j,i) = true;
            else
                k = any(PopObj(i,:)<PopObj(j,:)) - ...
                    any(PopObj(i,:)>PopObj(j,:));
                if k == 1
                    Dominate(i,j) = true;
                elseif k == -1
                    Dominate(j,i) = true;
                end
            end
        end
    end

    S = sum(Dominate,2);
    R = zeros(1,N);
    for i = 1:N
        R(i) = sum(S(Dominate(:,i)));
    end

    if N == 1
        D = 0;
    else
        Geometry = opfpim3e_scale_geometry(PopObj);
        Distance = pdist2(Geometry,Geometry);
        Distance(logical(eye(N))) = inf;
        Distance = sort(Distance,2);
        kth = max(1,min(N-1,floor(sqrt(N))));
        D = 1./(Distance(:,kth)+2);
    end

    Fitness = R + D';
end





