function [Repaired,diag] = opfpim3e_legacy_branch_repair(...
    Problem,Candidates)
%OPFPI_C3_BRANCH_AWARE_REPAIR Two-path branch-aware state restoration.
%
% For each selected infeasible candidate:
%   Path 1 starts from the inherited evolved state.
%   Path 2 starts from a deterministic flat-voltage state while preserving
%          the exact same controls.
%
% Each path uses the frozen C2 signed residual, analytic state Jacobian,
% spectral damping, trust scaling, and canonical halving backtracking.
%
% Selection:
%   - feasible beats infeasible;
%   - among infeasible candidates, lower physical residual key wins;
%   - among strictly feasible candidates, lower VD then higher min|V| wins.
%
% All calls to Problem.Evaluation, including flat-start evaluations, count
% in the real FE budget. Analytic residual/Jacobian calls are disclosed
% separately and do not consume FE.

    diag = initialize_diag();
    Repaired = [];

    if isempty(Candidates) || Problem.FE >= Problem.maxFE
        return;
    end

    tol = 1e-4;
    stateDim = 26;
    maxNewton = max(1,ceil(log10(1/tol)));
    maxBacktrack = maxNewton;
    maxCandidates = min(ceil(sqrt(Problem.N)),length(Candidates));

    infeas = find(~all(Candidates.cons<=0,2));
    if isempty(infeas)
        return;
    end

    DecInfeas = Candidates(infeas).decs;
    Hall = opfpim3e_legacy_signed_residual_jacobian(DecInfeas,false);
    diag.legacyResidualCalls = diag.legacyResidualCalls + size(Hall,1);

    absH = abs(Hall);
    profile = [max(absH,[],2),sqrt(sum(Hall.^2,2)),sum(absH,2)];
    score = opfpim3e_rank_sum(opfpim3e_normalize_columns(profile));
    [~,ord] = sort(score,'ascend');
    chosen = infeas(ord(1:min(maxCandidates,numel(ord))));

    % Keep the same per-generation maximum repair FE envelope as C2.
    remainingAtEntry = max(0,Problem.maxFE-Problem.FE);
    repairEvalEnvelope = min(remainingAtEntry,...
        (stateDim+maxBacktrack)*maxCandidates);

    stateRange = Problem.upper(1:stateDim)-Problem.lower(1:stateDim);
    evalStart = Problem.FE;

    for ii = 1:numel(chosen)
        if exhausted(Problem,evalStart,repairEvalEnvelope)
            break;
        end

        base = Candidates(chosen(ii));
        diag.repairCalls = diag.repairCalls + 1;
        diag.inheritedStarts = diag.inheritedStarts + 1;
        if isnan(diag.firstRepairCallFE)
            diag.firstRepairCallFE = Problem.FE;
        end

        [inherited,pi] = run_path(Problem,base,stateRange,...
            maxNewton,maxBacktrack,evalStart,repairEvalEnvelope);
        diag = merge_path_diag(diag,pi);

        flat = [];
        if ~exhausted(Problem,evalStart,repairEvalEnvelope)
            flatDec = opfpim3e_legacy_flat_start_dec(base.dec,Problem);
            flat = Problem.Evaluation(flatDec);
            diag.trialEvals = diag.trialEvals + 1;
            diag.evals = diag.evals + 1;
            diag.referenceStarts = diag.referenceStarts + 1;
            diag.referenceStartEvals = diag.referenceStartEvals + 1;

            if ~all(flat.cons<=0) && ...
                    ~exhausted(Problem,evalStart,repairEvalEnvelope)
                [flat,pf] = run_path(Problem,flat,stateRange,...
                    maxNewton,maxBacktrack,evalStart,repairEvalEnvelope);
                diag = merge_path_diag(diag,pf);
            end
        end

        [selected,source,usedBranchTie,extraResidualCalls] = ...
            choose_path(inherited,flat);
        diag.legacyResidualCalls = ...
            diag.legacyResidualCalls + extraResidualCalls;

        if usedBranchTie
            diag.branchTieBreakWins = diag.branchTieBreakWins + 1;
        end

        [betterThanBase,extraResidualCalls] = ...
            is_better_than_base(selected,base);
        diag.legacyResidualCalls = ...
            diag.legacyResidualCalls + extraResidualCalls;

        if betterThanBase
            Repaired = [Repaired,selected]; %#ok<AGROW>
            diag.acceptedCandidates = diag.acceptedCandidates + 1;

            if source == 1
                diag.inheritedPathWins = diag.inheritedPathWins + 1;
                if all(selected.cons<=0)
                    diag.strictFeasibleFromInherited = ...
                        diag.strictFeasibleFromInherited + 1;
                end
            elseif source == 2
                diag.referencePathWins = diag.referencePathWins + 1;
                if all(selected.cons<=0)
                    diag.strictFeasibleFromFlat = ...
                        diag.strictFeasibleFromFlat + 1;
                end
            end

            if all(selected.cons<=0)
                bk = opfpim3e_legacy_branch_key(selected.dec);
                diag.acceptedFeasibleVoltageDeviationSum = ...
                    diag.acceptedFeasibleVoltageDeviationSum + bk(1);
                diag.acceptedFeasibleMinVoltageSum = ...
                    diag.acceptedFeasibleMinVoltageSum - bk(2);
                diag.acceptedFeasibleCount = ...
                    diag.acceptedFeasibleCount + 1;
            end
        end
    end
end

function [current,d] = run_path(Problem,current,stateRange,...
    maxNewton,maxBacktrack,evalStart,repairEvalEnvelope)

    d = initialize_path_diag();
    currentDec = current.dec;
    h = opfpim3e_legacy_signed_residual_jacobian(currentDec,false);
    d.legacyResidualCalls = d.legacyResidualCalls + 1;
    currentKey = residual_key(h);

    for it = 1:maxNewton
        if all(current.cons<=0) || ...
                exhausted(Problem,evalStart,repairEvalEnvelope)
            break;
        end

        [h,J] = opfpim3e_legacy_signed_residual_jacobian(currentDec,true);
        d.legacyResidualCalls = d.legacyResidualCalls + 1;
        d.legacyAnalyticJacobianCalls = d.legacyAnalyticJacobianCalls + 1;
        d.newtonIterations = d.newtonIterations + 1;

        J = J(:,:,1);
        if any(~isfinite(h(:))) || any(~isfinite(J(:)))
            break;
        end

        [U,S,V] = svd(J,'econ');
        s = builtin('diag',S);
        if isempty(s) || max(s)<=eps
            break;
        end

        smax = max(s);
        smin = min(s);
        if smin<=eps*smax
            condJ = Inf;
        else
            condJ = smax/smin;
        end

        d.conditionNumberCount = d.conditionNumberCount + 1;
        if isfinite(condJ)
            d.conditionNumberSum = d.conditionNumberSum + condJ;
            d.conditionNumberMax = max(d.conditionNumberMax,condJ);
        else
            d.conditionNumberMax = Inf;
        end

        kappaTarget = eps^(-1/2);
        lambda = max(0,...
            (smax^2-kappaTarget*smin^2)/(kappaTarget-1));
        filt = s./(s.^2+lambda);
        delta = -V*(filt.*(U'*h'));

        if any(~isfinite(delta))
            break;
        end

        rel = max(abs(delta')./max(stateRange,eps));
        if rel>1
            delta = delta/rel;
        end

        acceptedThisIteration = false;
        alphas = 2.^-(0:maxBacktrack-1);

        for a = alphas
            if exhausted(Problem,evalStart,repairEvalEnvelope)
                break;
            end

            trialDec = currentDec;
            trialDec(1:26) = trialDec(1:26)+a*delta';
            trialDec = max(min(trialDec,Problem.upper),Problem.lower);

            trial = Problem.Evaluation(trialDec);
            d.trialEvals = d.trialEvals + 1;
            d.evals = d.evals + 1;
            d.lineSearchTrials = d.lineSearchTrials + 1;

            ht = opfpim3e_legacy_signed_residual_jacobian(trialDec,false);
            d.legacyResidualCalls = d.legacyResidualCalls + 1;
            trialKey = residual_key(ht);

            if opfpim3e_lexicographic_less(trialKey,currentKey)
                current = trial;
                currentDec = trialDec;
                currentKey = trialKey;
                acceptedThisIteration = true;
                d.acceptedSteps = d.acceptedSteps + 1;
                break;
            end
        end

        if ~acceptedThisIteration
            break;
        end
    end
end

function [selected,source,branchTie,residualCalls] = ...
    choose_path(inherited,flat)

    selected = inherited;
    source = 1;
    branchTie = false;
    residualCalls = 0;

    if isempty(flat)
        return;
    end

    fi = all(inherited.cons<=0);
    ff = all(flat.cons<=0);

    if ff && ~fi
        selected = flat;
        source = 2;
        return;
    elseif fi && ~ff
        return;
    elseif fi && ff
        ki = opfpim3e_legacy_branch_key(inherited.dec);
        kf = opfpim3e_legacy_branch_key(flat.dec);

        if opfpim3e_lexicographic_less(kf,ki)
            selected = flat;
            source = 2;
            branchTie = true;
        end
        return;
    end

    hi = opfpim3e_legacy_signed_residual_jacobian(inherited.dec,false);
    hf = opfpim3e_legacy_signed_residual_jacobian(flat.dec,false);
    residualCalls = residualCalls + 2;
    ki = residual_key(hi);
    kf = residual_key(hf);

    if opfpim3e_lexicographic_less(kf,ki)
        selected = flat;
        source = 2;
    elseif ~opfpim3e_lexicographic_less(ki,kf)
        bi = opfpim3e_legacy_branch_key(inherited.dec);
        bf = opfpim3e_legacy_branch_key(flat.dec);
        if opfpim3e_lexicographic_less(bf,bi)
            selected = flat;
            source = 2;
            branchTie = true;
        end
    end
end

function [tf,residualCalls] = is_better_than_base(candidate,base)
    residualCalls = 0;
    fc = all(candidate.cons<=0);
    fb = all(base.cons<=0);

    if fc && ~fb
        tf = true;
        return;
    elseif ~fc && fb
        tf = false;
        return;
    elseif fc && fb
        tf = opfpim3e_lexicographic_less(...
            opfpim3e_legacy_branch_key(candidate.dec),...
            opfpim3e_legacy_branch_key(base.dec));
        return;
    end

    hc = opfpim3e_legacy_signed_residual_jacobian(candidate.dec,false);
    hb = opfpim3e_legacy_signed_residual_jacobian(base.dec,false);
    residualCalls = residualCalls + 2;
    tf = opfpim3e_lexicographic_less(...
        residual_key(hc),residual_key(hb));
end

function key = residual_key(h)
    ah = abs(h(:));
    key = [max(ah),sqrt(sum(h(:).^2)),sum(ah)];
end

function tf = exhausted(Problem,evalStart,envelope)
    tf = Problem.FE>=Problem.maxFE || ...
        Problem.FE-evalStart>=envelope;
end

function d = initialize_diag()
    d = struct(...
        'repairCalls',0,...
        'acceptedCandidates',0,...
        'acceptedSteps',0,...
        'trialEvals',0,...
        'evals',0,...
        'legacyResidualCalls',0,...
        'legacyAnalyticJacobianCalls',0,...
        'newtonIterations',0,...
        'lineSearchTrials',0,...
        'conditionNumberSum',0,...
        'conditionNumberCount',0,...
        'conditionNumberMax',0,...
        'firstRepairCallFE',NaN,...
        'inheritedStarts',0,...
        'referenceStarts',0,...
        'referenceStartEvals',0,...
        'inheritedPathWins',0,...
        'referencePathWins',0,...
        'branchTieBreakWins',0,...
        'strictFeasibleFromInherited',0,...
        'strictFeasibleFromFlat',0,...
        'acceptedFeasibleVoltageDeviationSum',0,...
        'acceptedFeasibleMinVoltageSum',0,...
        'acceptedFeasibleCount',0);
end

function d = initialize_path_diag()
    d = struct(...
        'acceptedSteps',0,...
        'trialEvals',0,...
        'evals',0,...
        'legacyResidualCalls',0,...
        'legacyAnalyticJacobianCalls',0,...
        'newtonIterations',0,...
        'lineSearchTrials',0,...
        'conditionNumberSum',0,...
        'conditionNumberCount',0,...
        'conditionNumberMax',0);
end

function d = merge_path_diag(d,p)
    names = fieldnames(p);
    for i = 1:numel(names)
        name = names{i};
        if strcmp(name,'conditionNumberMax')
            d.(name) = max(d.(name),p.(name));
        else
            d.(name) = d.(name)+p.(name);
        end
    end
end





