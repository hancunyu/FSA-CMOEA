function Y = opfpim3e_scale_geometry(X)
%OPFPIU_SCALE_GEOMETRY Robust column-wise min-max scaling for geometry.
%
% This function is used only by distance-, density- and angle-based search
% components. Pareto dominance and the original objective values are not
% changed. Constant columns are mapped to zero. Non-finite values are
% guarded deterministically, although normal OPF evaluations should be
% finite before reaching selection.

    Y = zeros(size(X));
    if isempty(X)
        return;
    end

    for j = 1:size(X,2)
        column = X(:,j);
        finiteMask = isfinite(column);

        if ~any(finiteMask)
            continue;
        end

        lo = min(column(finiteMask));
        hi = max(column(finiteMask));
        scale = hi-lo;
        threshold = 64*eps*max([1,abs(lo),abs(hi)]);

        if scale <= threshold
            Y(finiteMask,j) = 0;
        else
            Y(finiteMask,j) = (column(finiteMask)-lo)./scale;
        end

        % Deterministic guards for unexpected non-finite values.
        Y(isinf(column) & column>0,j) = 1;
        Y(isinf(column) & column<0,j) = 0;
        Y(isnan(column),j) = 0.5;
    end
end





