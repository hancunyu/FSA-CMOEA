function diag = opfpim3e_diag_init(Algorithm,Problem,profile,P1,P2)
    diag = struct();
    diag.algorithm = class(Algorithm);
    diag.problem = class(Problem);
    diag.mode = char(profile.mode);
    diag.run = Algorithm.run;
    diag.initialPopulation1SHA = opfpim3e_sha256_matrix(P1.decs);
    diag.initialPopulation2SHA = opfpim3e_sha256_matrix(P2.decs);

    diag.generationCount = 0;
    diag.evolutionEvals = 0;
    diag.repairEvals = 0;
    diag.repairActiveGenerations = 0;
    diag.repairSkippedGenerations = 0;
    diag.repairGateFE = NaN;
    diag.relaxationProgressFinal = 0;

    diag.repairCalls = 0;
    diag.acceptedCandidates = 0;
    diag.acceptedSteps = 0;
    diag.trialEvals = 0;
    diag.referenceStarts = 0;
    diag.referenceStartEvals = 0;
    diag.inheritedPathWins = 0;
    diag.referencePathWins = 0;
    diag.branchTieBreakWins = 0;
    diag.learnedSensitivityCalls = 0;
    diag.localModelFailures = 0;
    diag.legacyResidualCalls = 0;
    diag.legacyAnalyticJacobianCalls = 0;
    diag.conditionNumberSum = 0;
    diag.conditionNumberCount = 0;
    diag.conditionNumberMax = 0;
    diag.firstRepairCallFE = NaN;
    diag.odpCalls = 0;
    diag.odpScoreSum = 0;
    diag.odpPenaltySum = 0;
    diag.odpMeanScore = NaN;
    diag.odpMeanPenalty = NaN;
    diag.odpDemotedCount = 0;
    diag.odpAcceptedCount = 0;
    diag.odpFeasibleAccepted = 0;
    diag.odpFeasibleDemoted = 0;
    diag.repairObjWorsenSum = 0;
    diag.repairObjWorsenMean = NaN;
    diag.repairObjWorsenMax = 0;

    diag.trace = struct(...
        'FE',[],...
        'generation',[],...
        'relaxationProgress',[],...
        'cumulativeEvolutionEvals',[],...
        'cumulativeRepairEvals',[],...
        'repairGateActive',[],...
        'feasible1',[],...
        'feasible2',[],...
        'publishedFeasible',[],...
        'bestCVsum',[],...
        'bestCVmax',[],...
        'population1SHA',{{}},...
        'population2SHA',{{}});
end





