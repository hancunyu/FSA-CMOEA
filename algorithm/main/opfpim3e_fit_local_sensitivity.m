function [delta,info] = opfpim3e_fit_local_sensitivity(...
    currentDec,currentCon,Model,active,Problem)
%OPFPIU_FIT_LOCAL_SENSITIVITY Estimate d(g_active)/d(x) from evaluated data.
%
% No new problem evaluation is performed. The local linear model uses the
% nearest already evaluated decisions and spectral damping derived from
% machine precision, not a tuned regularization constant.

    info = struct('success',false,'conditionNumber',NaN);
    delta = zeros(1,Problem.D);

    if isempty(active) || size(Model.Dec,1) < 3
        return;
    end

    range = Problem.upper-Problem.lower;
    range(range<=eps) = 1;

    X = (Model.Dec-Problem.lower)./range;
    x0 = (currentDec-Problem.lower)./range;

    dist = sqrt(sum((X-x0).^2,2));
    valid = isfinite(dist) & dist > 32*eps;

    nValid = sum(valid);
    if nValid < 2
        return;
    end

    kLocal = min(nValid,max(Problem.D+1,ceil(sqrt(nValid))));
    validIdx = find(valid);
    [~,order] = sort(dist(validIdx),'ascend');
    idx = validIdx(order(1:kLocal));

    dX = X(idx,:)-x0;
    G = Model.Con(:,active);
    g0 = currentCon(active);

    gLo = min(G,[],1);
    gHi = max(G,[],1);
    gScale = max([gHi-gLo;abs(g0);ones(1,numel(active))*eps],[],1);
    gScale(gScale<=eps) = 1;

    dG = (Model.Con(idx,active)-g0)./gScale;
    target = g0./gScale;

    finiteRows = all(isfinite(dX),2) & all(isfinite(dG),2);
    dX = dX(finiteRows,:);
    dG = dG(finiteRows,:);

    if size(dX,1) < 2 || norm(dX,'fro') <= eps
        return;
    end

    [B,condFit,okFit] = local_damped_least_squares(dX,dG);
    if ~okFit || any(~isfinite(B(:)))
        return;
    end

    A = B';
    [stepColumn,condSolve,okSolve] = ...
        local_damped_least_squares(A,-target');

    if ~okSolve || any(~isfinite(stepColumn(:)))
        return;
    end

    delta = stepColumn';
    rel = max(abs(delta));
    if rel > 1
        delta = delta/rel;
    end

    if norm(delta) <= 64*eps
        return;
    end

    info.success = true;
    info.conditionNumber = max(condFit,condSolve);
end

function [X,condA,ok] = local_damped_least_squares(A,B)
% Solve min ||A*X-B|| with spectral damping.

    X = zeros(size(A,2),size(B,2));
    condA = Inf;
    ok = false;

    if isempty(A) || any(~isfinite(A(:))) || any(~isfinite(B(:)))
        return;
    end

    [U,S,V] = svd(A,'econ');
    s = diag(S);

    if isempty(s) || max(s) <= eps
        return;
    end

    smax = max(s);
    smin = min(s);

    if smin <= eps*smax
        condA = Inf;
    else
        condA = smax/smin;
    end

    kappaTarget = eps^(-1/2);
    lambda = max(0,...
        (smax^2-kappaTarget*smin^2)/(kappaTarget-1));

    filter = s./(s.^2+lambda);
    X = V*(filter.*(U'*B));
    ok = all(isfinite(X(:)));
end





