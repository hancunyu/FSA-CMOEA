function [RepairedMain,RepairedAux,diag] = opfpim3e_control_branch_repair(...
    Problem,Samples,Candidates,profile)
%opfpim3e_CONTROL_BRANCH_REPAIR Control-space residual restoration with ODP.
%
% The evolved candidate is always the primary path. A deterministic
% reference-control path is used only as a fallback when the inherited path
% remains infeasible. This differs from the state-encoded RWMOP case, where
% the two C3 branches share identical control variables.
%
% Selection rules:
%   1) feasible beats infeasible;
%   2) among infeasible paths, lower constraint-violation key wins;
%   3) when both paths are feasible, the path with the smaller normalized
%      decision disturbance from the evolved candidate wins;
%   4) no objective, including LVD, is used to choose a repair branch.
%
% Every evaluated trial passes through Problem.Evaluation and consumes FE.

    diag = local_init_diag();
    RepairedMain = [];
    RepairedAux = [];

    if isempty(Candidates) || isempty(Samples) || ...
            Problem.FE >= Problem.maxFE
        return;
    end

    Model.Dec = Samples.decs;
    Model.Con = Samples.cons;
    Model.Obj = Samples.objs;

    finiteRows = all(isfinite(Model.Dec),2) & ...
        all(isfinite(Model.Con),2) & all(isfinite(Model.Obj),2);
    Model.Dec = Model.Dec(finiteRows,:);
    Model.Con = Model.Con(finiteRows,:);
    Model.Obj = Model.Obj(finiteRows,:);

    if isempty(Model.Dec)
        return;
    end

    [~,uniqueIdx] = unique(Model.Dec,'rows','stable');
    Model.Dec = Model.Dec(uniqueIdx,:);
    Model.Con = Model.Con(uniqueIdx,:);
    Model.Obj = Model.Obj(uniqueIdx,:);

    objLo = min(Model.Obj,[],1);
    objScale = max(Model.Obj,[],1)-objLo;
    objScale(objScale<=eps) = 1;

    infeasible = find(~all(Candidates.cons<=0,2));
    if isempty(infeasible)
        return;
    end

    features = opfpim3e_residual_features(...
        Candidates(infeasible).cons,profile);
    score = opfpim3e_rank_sum(opfpim3e_normalize_columns(features));
    [~,order] = sort(score,'ascend');

    maxCandidates = min(ceil(sqrt(Problem.N)),numel(infeasible));
    chosen = infeasible(order(1:maxCandidates));

    tol = 1e-4;
    maxNewton = max(1,ceil(log10(1/tol)));
    maxBacktrack = maxNewton;

    remainingAtEntry = max(0,Problem.maxFE-Problem.FE);

    % In the control-only representation, restoration must not consume more
    % than one population-equivalent FE block in a generation. This is a
    % deterministic structural limit, not a tuned problem parameter.
    repairEnvelope = min([...
        remainingAtEntry,...
        Problem.N,...
        (Problem.D+maxBacktrack)*maxCandidates]);

    evalStart = Problem.FE;

    for i = 1:numel(chosen)
        if local_exhausted(Problem,evalStart,repairEnvelope)
            break;
        end

        base = Candidates(chosen(i));
        diag.repairCalls = diag.repairCalls+1;
        if isnan(diag.firstRepairCallFE)
            diag.firstRepairCallFE = Problem.FE;
        end

        % Primary path: preserve the evolved control vector and make the
        % smallest residual-guided control correction available.
        [inherited,Model,pathDiag] = local_run_path(...
            Problem,base,Model,profile,objLo,objScale,...
            maxNewton,maxBacktrack,evalStart,repairEnvelope);
        diag = local_merge(diag,pathDiag);

        reference = [];

        % The reference-control path is a fallback only. If the inherited
        % path is already feasible, evaluating a second control vector would
        % spend FE and may distort the evolved objective-search direction.
        if ~all(inherited.cons<=0,2) && ...
                ~local_exhausted(Problem,evalStart,repairEnvelope)

            referenceDec = opfpim3e_reference_start_dec(...
                base.dec,Problem,profile);

            if any(abs(referenceDec-base.dec) > 64*eps)
                reference = Problem.Evaluation(referenceDec);
                diag.trialEvals = diag.trialEvals+1;
                diag.referenceStarts = diag.referenceStarts+1;
                diag.referenceStartEvals = ...
                    diag.referenceStartEvals+1;
                Model = local_append_model(Model,reference);

                if ~all(reference.cons<=0,2) && ...
                        ~local_exhausted(Problem,evalStart,repairEnvelope)
                    [reference,Model,pathDiag] = local_run_path(...
                        Problem,reference,Model,profile,objLo,objScale,...
                        maxNewton,maxBacktrack,evalStart,repairEnvelope);
                    diag = local_merge(diag,pathDiag);
                end
            end
        end

        [selected,source,usedDistanceTieBreak] = ...
            local_choose_path(inherited,reference,base,...
                Problem,profile,objLo,objScale);

        if usedDistanceTieBreak
            diag.branchTieBreakWins = diag.branchTieBreakWins+1;
        end

        if opfpim3e_solution_better(selected,base,...
                profile,objLo,objScale)

            [directionPreserved,odpDiag] = local_odp_direction_score(...
                Problem,profile,base,selected,Model,objLo,objScale);
            diag = local_merge(diag,odpDiag);

            % M3E: repair-to-auxiliary exposure is value-aware. Repairs
            % that preserve the main direction are always exposed; demoted
            % repairs enter the relaxed auxiliary pool only when they carry
            % sparse-frontier or engineering-extreme value. This keeps the
            % auxiliary search from being over-dominated by generic repair
            % traces while still protecting useful OPF extreme points.
            diag.acceptedCandidates = diag.acceptedCandidates+1;

            if directionPreserved
                RepairedMain = [RepairedMain,selected]; %#ok<AGROW>
                RepairedAux = [RepairedAux,selected]; %#ok<AGROW>
                diag.odpAcceptedCount = diag.odpAcceptedCount+1;
                if all(selected.cons<=0,2)
                    diag.odpFeasibleAccepted = diag.odpFeasibleAccepted+1;
                end
            else
                if local_aux_value_worthy(selected,base,Model,objLo,objScale)
                    RepairedAux = [RepairedAux,selected]; %#ok<AGROW>
                end
                diag.odpDemotedCount = diag.odpDemotedCount+1;
                if all(selected.cons<=0,2)
                    diag.odpFeasibleDemoted = diag.odpFeasibleDemoted+1;
                end
            end

            if source == 1
                diag.inheritedPathWins = diag.inheritedPathWins+1;
            else
                diag.referencePathWins = diag.referencePathWins+1;
            end
        end
    end
end

function [selected,source,usedDistanceTieBreak] = local_choose_path(...
    inherited,reference,base,Problem,profile,objLo,objScale)
% Choose the repair path without introducing objective preference.

    selected = inherited;
    source = 1;
    usedDistanceTieBreak = false;

    if isempty(reference)
        return;
    end

    inheritedFeasible = all(inherited.cons<=0,2);
    referenceFeasible = all(reference.cons<=0,2);

    if referenceFeasible && ~inheritedFeasible
        selected = reference;
        source = 2;
        return;
    elseif inheritedFeasible && ~referenceFeasible
        return;
    elseif ~inheritedFeasible && ~referenceFeasible
        if opfpim3e_solution_better(reference,inherited,...
                profile,objLo,objScale)
            selected = reference;
            source = 2;
        end
        return;
    end

    % Both paths are feasible. Select the least disruptive correction in
    % normalized decision space. Exact ties retain the inherited path.
    range = Problem.upper-Problem.lower;
    range(range<=eps) = 1;

    inheritedDistance = norm(...
        (inherited.dec-base.dec)./range,2);
    referenceDistance = norm(...
        (reference.dec-base.dec)./range,2);

    distanceTolerance = 64*eps*max(...
        [1,inheritedDistance,referenceDistance]);

    if referenceDistance < inheritedDistance-distanceTolerance
        selected = reference;
        source = 2;
        usedDistanceTieBreak = true;
    end
end

function [current,Model,diag] = local_run_path(...
    Problem,current,Model,profile,objLo,objScale,...
    maxNewton,maxBacktrack,evalStart,repairEnvelope)

    diag = local_init_path_diag();

    for iteration = 1:maxNewton %#ok<NASGU>
        if all(current.cons<=0,2) || ...
                local_exhausted(Problem,evalStart,repairEnvelope)
            break;
        end

        active = find(current.cons>0);
        [delta,modelInfo] = opfpim3e_fit_local_sensitivity(...
            current.dec,current.cons,Model,active,Problem);
        diag.learnedSensitivityCalls = ...
            diag.learnedSensitivityCalls+1;

        if isfinite(modelInfo.conditionNumber)
            diag.conditionNumberCount = ...
                diag.conditionNumberCount+1;
            diag.conditionNumberSum = ...
                diag.conditionNumberSum+modelInfo.conditionNumber;
            diag.conditionNumberMax = max(...
                diag.conditionNumberMax,modelInfo.conditionNumber);
        elseif isinf(modelInfo.conditionNumber)
            diag.conditionNumberMax = Inf;
        end

        if ~modelInfo.success
            diag.localModelFailures = diag.localModelFailures+1;
            delta = local_fallback_direction(current,Model,Problem);
        end

        if norm(delta) <= 64*eps || any(~isfinite(delta))
            break;
        end

        accepted = false;
        range = Problem.upper-Problem.lower;
        range(range<=eps) = 1;
        x0 = (current.dec-Problem.lower)./range;

        for alpha = 2.^-(0:maxBacktrack-1)
            if local_exhausted(Problem,evalStart,repairEnvelope)
                break;
            end

            trialNorm = x0+alpha*delta;
            trialDec = Problem.lower+trialNorm.*range;
            trialDec = opfpim3e_project_dec(...
                Problem,trialDec,profile);

            if all(abs(trialDec-current.dec) <= 64*eps)
                continue;
            end

            trial = Problem.Evaluation(trialDec);
            diag.trialEvals = diag.trialEvals+1;
            Model = local_append_model(Model,trial);

            if opfpim3e_solution_better(trial,current,...
                    profile,objLo,objScale)
                current = trial;
                accepted = true;
                diag.acceptedSteps = diag.acceptedSteps+1;
                break;
            end
        end

        if ~accepted
            break;
        end
    end
end

function [directionPreserved,diag] = local_odp_direction_score(Problem,profile,base,selected,Model,objLo,objScale)
%LOCAL_ODP_DIRECTION_SCORE Conservative balanced direction preservation.
% ODPB2 treats the original ODP rule as the baseline: if total normalized
% improvement already offsets total worsening, the repair is allowed into the
% main pool. Weak-objective and sparse-frontier signals only rescue repairs
% rejected by the baseline rule; they no longer veto baseline-acceptable
% trade-offs. This keeps the mechanism parameter-free and avoids over-steering
% the IEEE30_C4 front while still rewarding useful coverage expansion.

    diag = struct(...
        'odpCalls',1,...
        'odpScoreSum',0,...
        'odpPenaltySum',0,...
        'odpDemotedCount',0,...
        'odpAcceptedCount',0,...
        'odpFeasibleAccepted',0,...
        'odpFeasibleDemoted',0,...
        'repairObjWorsenSum',0,...
        'repairObjWorsenMax',0);

    if ~local_odp_guard_active(Problem,profile)
        directionPreserved = true;
        return;
    end

    scale = objScale;
    scale(scale <= eps) = 1;
    deltaObj = (selected.obj - base.obj) ./ scale;
    deltaObj(~isfinite(deltaObj)) = 0;

    worsening = max(0,deltaObj);
    improvement = max(0,-deltaObj);
    worsenMass = sum(worsening);
    improveMass = sum(improvement);

    modelObj = Model.Obj;
    modelObj = modelObj(all(isfinite(modelObj),2),:);
    if isempty(modelObj) || size(modelObj,2) ~= numel(deltaObj)
        modelObj = [base.obj;selected.obj];
    end

    z = (base.obj - objLo) ./ scale;
    z(~isfinite(z)) = 0;
    z = max(0,z);
    if max(z) > 0
        weakWeight = 1 + z ./ max(z);
    else
        weakWeight = ones(size(z));
    end

    weightedWorsen = sum(worsening .* weakWeight);
    weightedImprove = sum(improvement .* weakWeight);

    idealGain = sum(max(0,(objLo - selected.obj) ./ scale));
    sparseGain = local_frontier_sparse_gain(selected.obj,modelObj,objLo,scale);
    frontierCredit = idealGain + sparseGain;

    penalty = mean(worsening);
    score = (weightedImprove + frontierCredit) / ...
        (1 + weightedWorsen + weightedImprove + frontierCredit);

    diag.odpScoreSum = score;
    diag.odpPenaltySum = penalty;
    diag.repairObjWorsenSum = penalty;
    diag.repairObjWorsenMax = max(worsening);

    % Strict feasibility still has priority. A repaired candidate that is
    % already strictly feasible is allowed to enter the main pool, because
    % excluding feasible repairs can remove useful Pareto-extension points
    % and hurt HV/spacing in high-dimensional IEEE cases. The UNIFIED gate only
    % demotes accepted repairs that still carry feasibility risk.
    baselineTolerance = 64*eps*max(1,max(abs([worsenMass,improveMass])));
    balancedTolerance = 64*eps*max(1,max(abs([weightedWorsen,weightedImprove,frontierCredit])));
    baselinePass = worsenMass <= improveMass + baselineTolerance;
    balancedRescue = weightedWorsen <= weightedImprove + frontierCredit + balancedTolerance;
    feasibleRepair = all(selected.cons<=0,2);
    directionPreserved = feasibleRepair || baselinePass || balancedRescue;
end

function tf = local_odp_guard_active(Problem,profile)
% Activate objective-direction demotion only when the OPF control space is
% large enough for repair trials to dominate the search cost. Smaller OPF
% cases still record ODP diagnostics but keep every accepted repair in the
% strict main pool, preserving the original OPFPI behaviour inside the same
% algorithm.
    tf = ismember(profile.mode,["ieee30_control","ieee57_control"]) && ...
        Problem.D > ceil(sqrt(Problem.N))*3;
end

function tf = local_aux_value_worthy(selected,base,Model,objLo,objScale)
% Keep demoted repairs in the relaxed auxiliary pool only if they add
% identifiable Pareto-front value or improve at least one engineering
% objective beyond the current evaluated model envelope.

    scale = objScale;
    scale(scale <= eps) = 1;
    modelObj = Model.Obj;
    modelObj = modelObj(all(isfinite(modelObj),2),:);
    if isempty(modelObj) || size(modelObj,2) ~= numel(selected.obj)
        modelObj = [base.obj;selected.obj];
    end

    improvement = (base.obj - selected.obj) ./ scale;
    improvement(~isfinite(improvement)) = 0;
    singleObjectiveGain = max(improvement);

    modelLo = min(modelObj,[],1);
    extremeGain = max((modelLo - selected.obj) ./ scale);
    if ~isfinite(extremeGain)
        extremeGain = 0;
    end

    sparseGain = local_frontier_sparse_gain(selected.obj,modelObj,objLo,scale);
    tol = 64*eps*max(1,max(abs([singleObjectiveGain,extremeGain,sparseGain])));
    tf = singleObjectiveGain > tol || extremeGain > tol || sparseGain > tol;
end

function sparseGain = local_frontier_sparse_gain(obj,modelObj,objLo,objScale)
% Reward repairs that move into a sparse normalized objective region.

    sparseGain = 0;
    if size(modelObj,1) < 3
        return;
    end

    Z = (modelObj - objLo) ./ objScale;
    z = (obj - objLo) ./ objScale;
    Z(~isfinite(Z)) = 0;
    z(~isfinite(z)) = 0;

    d = sqrt(sum((Z - z).^2,2));
    d = d(isfinite(d) & d > 64*eps);
    if isempty(d)
        return;
    end

    nearest = min(d);
    scaleDist = median(d);
    if ~isfinite(scaleDist) || scaleDist <= 64*eps
        return;
    end

    sparseGain = max(0,nearest/scaleDist - 1);
end

function delta = local_fallback_direction(current,Model,Problem)
% Use only an already evaluated lower-violation neighbor. If none exists,
% return zero instead of pulling the inherited path toward a fixed control
% reference. The separate reference path remains available as fallback.

    range = Problem.upper-Problem.lower;
    range(range<=eps) = 1;
    X = (Model.Dec-Problem.lower)./range;
    x0 = (current.dec-Problem.lower)./range;

    V = max(0,Model.Con);
    sampleKey = [max(V,[],2),sum(V,2)];
    vc = max(0,current.cons);
    currentKey = [max(vc,[],2),sum(vc,2)];

    better = sampleKey(:,1) < currentKey(1)-64*eps | ...
        (abs(sampleKey(:,1)-currentKey(1))<=64*eps & ...
         sampleKey(:,2) < currentKey(2)-64*eps);

    if any(better)
        idx = find(better);
        distance = sqrt(sum((X(idx,:)-x0).^2,2));
        [~,best] = min(distance);
        delta = X(idx(best),:)-x0;
    else
        delta = zeros(1,Problem.D);
    end

    rel = max(abs(delta));
    if rel > 1
        delta = delta/rel;
    end
end

function Model = local_append_model(Model,solution)
    Model.Dec = [Model.Dec;solution.decs];
    Model.Con = [Model.Con;solution.cons];
    Model.Obj = [Model.Obj;solution.objs];
end

function tf = local_exhausted(Problem,evalStart,envelope)
    tf = Problem.FE>=Problem.maxFE || ...
        Problem.FE-evalStart>=envelope;
end

function diag = local_init_diag()
    diag = struct(...
        'repairCalls',0,...
        'acceptedCandidates',0,...
        'acceptedSteps',0,...
        'trialEvals',0,...
        'referenceStarts',0,...
        'referenceStartEvals',0,...
        'inheritedPathWins',0,...
        'referencePathWins',0,...
        'branchTieBreakWins',0,...
        'learnedSensitivityCalls',0,...
        'localModelFailures',0,...
        'conditionNumberSum',0,...
        'conditionNumberCount',0,...
        'conditionNumberMax',0,...
        'odpCalls',0,...
        'odpScoreSum',0,...
        'odpPenaltySum',0,...
        'odpMeanScore',NaN,...
        'odpMeanPenalty',NaN,...
        'odpDemotedCount',0,...
        'odpAcceptedCount',0,...
        'odpFeasibleAccepted',0,...
        'odpFeasibleDemoted',0,...
        'repairObjWorsenSum',0,...
        'repairObjWorsenMean',NaN,...
        'repairObjWorsenMax',0,...
        'firstRepairCallFE',NaN);
end

function diag = local_init_path_diag()
    diag = struct(...
        'acceptedSteps',0,...
        'trialEvals',0,...
        'learnedSensitivityCalls',0,...
        'localModelFailures',0,...
        'conditionNumberSum',0,...
        'conditionNumberCount',0,...
        'conditionNumberMax',0);
end

function total = local_merge(total,part)
    names = fieldnames(part);
    for k = 1:numel(names)
        name = names{k};
        if strcmp(name,'conditionNumberMax') || strcmp(name,'repairObjWorsenMax')
            total.(name) = max(total.(name),part.(name));
        else
            total.(name) = total.(name)+part.(name);
        end
    end
end





