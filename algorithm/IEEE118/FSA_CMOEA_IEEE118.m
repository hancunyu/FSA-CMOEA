classdef FSA_CMOEA_IEEE118 < ALGORITHM
    properties
        diagnostics = struct();
    end
% <2026> <multi> <real> <constrained>
% OPF physics-informed constrained multi-objective evolutionary algorithm.
% ODP: objective-direction-preserving physical restoration for IEEE30/57.
%
% Unified paper implementation.
%
% Evolutionary core:
%   1) two interacting populations;
%   2) DE/rand/1 and current-to-pbest/1 variation;
%   3) feasibility-first main-population selection;
%   4) relaxed-constraint auxiliary-population selection;
%   5) residual-guided physical restoration.
%
% Representation-specific restoration:
%   - RWMOP40--46 (D=34): frozen C3 signed-Jacobian restoration;
%   - IEEE30_MOOPF (D=24) and IEEE57_MOOPF (D=33): control-space local-sensitivity restoration.
%
% IEEE30 semantic corrections in this revision:
%   A) the auxiliary relaxation schedule advances by evolutionary
%      generations, as in the frozen C3 prototype, rather than by total FE;
%   B) once the published main population is entirely feasible, additional
%      control-space restoration is disabled and the remaining FE budget is
%      devoted to Pareto-front search;
%   C) feasible repair paths are compared by normalized disturbance from
%      the evolved candidate, not by any individual objective such as LVD;
%   D) all restored trials still pass through Problem.Evaluation and consume
%      real FE. No hidden evaluator or MATPOWER call is made by the algorithm.
%   E) all objective-space distances and angles are computed after online
%      min-max scaling, while Pareto dominance is evaluated on the original
%      objectives. This makes the search geometry invariant to engineering
%      units such as cost, loss, emission and voltage deviation.
%   F) the auxiliary relaxation progress is driven by the number of actual
%      evolutionary evaluations after the two initial populations, rather
%      than by a nominal one-population generation count.

    methods
        function main(Algorithm,Problem)
            profile = opfpim3e_validate_problem(Problem);

            Population1 = Problem.Initialization();
            Fitness1 = opfpim3e_calfitness(Population1.objs,Population1.cons);
            Zmin1 = min(Population1.objs,[],1);

            Population2 = Problem.Initialization();
            Fitness2 = opfpim3e_calfitness(Population2.objs,Population2.cons);
            Zmin2 = min(Population2.objs,[],1);

            cons0 = [Population1.cons;Population2.cons];
            cons0(cons0<0) = 0;
            VAR0 = max(sum(cons0,2));
            if ~isfinite(VAR0) || VAR0 <= 0
                VAR0 = 1;
            end

            % Use evolutionary-evaluation progress. Two initial populations
            % consume 2*N evaluations; repair evaluations do not accelerate
            % the auxiliary relaxation schedule.
            progress = 0;
            evolutionBudget = max(1,Problem.maxFE-2*Problem.N);
            generation = 0;

            diag = opfpim3e_diag_init(Algorithm,Problem,profile,...
                Population1,Population2);
            diag.revision = ...
                'unified_objective_direction_preserving_restoration';

            while true
                Published = Population1;
                diag = opfpim3e_diag_snapshot(diag,Problem,Population1,...
                    Population2,Published);

                traceIndex = numel(diag.trace.FE);
                diag.trace.generation(traceIndex,1) = generation;
                diag.trace.relaxationProgress(traceIndex,1) = progress;
                diag.trace.cumulativeEvolutionEvals(traceIndex,1) = ...
                    diag.evolutionEvals;
                diag.trace.cumulativeRepairEvals(traceIndex,1) = ...
                    diag.repairEvals;
                diag.trace.repairGateActive(traceIndex,1) = ...
                    double(opfpim3e_repair_enabled(Population1,profile));

                opfpim3e_save_diag(Algorithm,Problem,diag);

                if ~Algorithm.NotTerminated(Published)
                    break;
                end

                remaining = max(0,Problem.maxFE-Problem.FE);
                if remaining <= 0
                    break;
                end

                cp = (-log(VAR0)-6)/log(1-0.5);
                VAR = VAR0*max(0,1-progress)^cp;

                % Keep the frozen C3 two-population offspring allocation.
                nOff1 = min(Problem.N,ceil(remaining/2));
                nOff2 = min(Problem.N,floor(remaining/2));

                generation = generation+1;
                diag.generationCount = generation;

                evolutionFEStart = Problem.FE;

                Offspring1 = opfpim3e_generate_offspring(Problem,...
                    Population1,Fitness1,Zmin1,nOff1,profile);
                Offspring2 = opfpim3e_generate_offspring(Problem,...
                    Population2,Fitness2,Zmin2,nOff2,profile);

                diag.evolutionEvals = diag.evolutionEvals + ...
                    (Problem.FE-evolutionFEStart);

                if ~isempty(Offspring1)
                    Zmin1 = min([Zmin1;Offspring1.objs],[],1);
                end
                if ~isempty(Offspring2)
                    Zmin2 = min([Zmin2;Offspring2.objs],[],1);
                end

                Samples = [Population1,Population2,Offspring1,Offspring2];
                Candidates = [Offspring1,Offspring2];
                RepairedMain = [];
                RepairedAux = [];

                repairEnabled = opfpim3e_repair_enabled(Population1,profile);
                repairFEStart = Problem.FE;

                if repairEnabled && Problem.FE < Problem.maxFE
                    diag.repairActiveGenerations = ...
                        diag.repairActiveGenerations+1;

                    [RepairedMain,RepairedAux,repairDiag] = opfpim3e_branch_aware_repair(...
                        Problem,Samples,Candidates,profile);
                    diag = opfpim3e_merge_diag(diag,repairDiag);
                    if isfield(diag,'odpCalls') && diag.odpCalls > 0
                        if isfield(diag,'odpScoreSum')
                            diag.odpMeanScore = ...
                                diag.odpScoreSum/max(1,diag.odpCalls);
                        end
                        if isfield(diag,'odpPenaltySum')
                            diag.odpMeanPenalty = ...
                                diag.odpPenaltySum/max(1,diag.odpCalls);
                        end
                        if isfield(diag,'repairObjWorsenSum')
                            diag.repairObjWorsenMean = ...
                                diag.repairObjWorsenSum/max(1,diag.odpCalls);
                        end
                    end
                else
                    diag.repairSkippedGenerations = ...
                        diag.repairSkippedGenerations+1;

                    if ismember(profile.mode,["ieee30_control","ieee57_control","ieee118_control"]) && ...
                            isnan(diag.repairGateFE)
                        diag.repairGateFE = Problem.FE;
                    end
                end

                diag.repairEvals = diag.repairEvals + ...
                    (Problem.FE-repairFEStart);

                % ODP demotes objective-damaging repaired candidates from
                % the strict main pool, while still exposing them to the
                % relaxed auxiliary pool. This changes no objective value and
                % adds no extra power-flow evaluation.
                Pool1 = [Population1,Offspring1,Offspring2,RepairedMain];
                Pool2 = [Population2,Offspring2,Offspring1,RepairedAux];

                [Population1,Fitness1] = ...
                    opfpim3e_main_selection(Pool1,Problem.N,true);
                [Population2,Fitness2] = ...
                    opfpim3e_aux_selection_original(Pool2,Problem.N,VAR);

                progress = min(1,diag.evolutionEvals/evolutionBudget);
                diag.relaxationProgressFinal = progress;
            end
        end
    end
end





