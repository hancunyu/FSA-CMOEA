function key = opfpim3e_solution_key(solution,profile,objLo,objScale)
%OPFPIU_SOLUTION_KEY Objective-neutral feasibility key for restoration.
% Lower is better.
%
% For infeasible solutions, the key is based only on constraint violation.
% For feasible solutions, an aggregate normalized objective key is retained
% only as a deterministic generic fallback. No individual objective (such
% as IEEE30 LVD) receives priority in the restoration operator.

    con = solution.cons;
    obj = solution.objs;
    feasible = all(con<=0,2);

    if feasible
        on = (obj-objLo)./objScale;
        on(~isfinite(on)) = 0;
        key = [0,sum(on),on];
    else
        V = max(0,con);
        base = [max(V,[],2),sum(V,2)];
        groups = profile.constraintGroups;
        grouped = [];
        for k = 1:numel(groups)
            idx = groups{k};
            idx = idx(idx>=1 & idx<=numel(V));
            if ~isempty(idx)
                grouped = [grouped,sum(V(idx)),max(V(idx))]; %#ok<AGROW>
            end
        end
        key = [1,base,grouped];
    end
end





