function result = run_fsa_case(caseName, options)
%RUN_FSA_CASE Run one FSA-CMOEA example with the paper's default settings.
% Set PLATEMO_ROOT and MATPOWER_ROOT as described in README.md.
arguments
    caseName (1,1) string
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 1
    options.MaxFE (1,1) double {mustBeInteger,mustBePositive} = 50000
end
    root = fileparts(mfilename('fullpath'));
    caseIndex = find(caseName == "Case" + string(1:8), 1);
    isIEEE118 = caseName == "IEEE118";
    isRWMOP = any(caseName == "RWMOP" + string(40:46));
    assert(~isempty(caseIndex) || isIEEE118 || isRWMOP, ...
        'FSA:UnknownCase', 'Use Case1--Case8, IEEE118, or RWMOP40--RWMOP46.');
    assert(options.MaxFE>200, 'FSA:Budget', ...
        'MaxFE must exceed 200 evaluations for the two initial populations.');
    originalPath = path;
    restorePath = onCleanup(@() path(originalPath)); %#ok<NASGU>
    platemoRoot = getenv('PLATEMO_ROOT');
    assert(isfile(fullfile(platemoRoot,'Algorithms','ALGORITHM.m')), ...
        'FSA:PlatEMO', 'Set PLATEMO_ROOT to the PlatEMO directory.');
    addpath(genpath(platemoRoot), '-end');
    addpath(fullfile(platemoRoot,'Algorithms'), fullfile(platemoRoot,'Problems'), '-begin');
    if ~isRWMOP
        matpowerRoot = getenv('MATPOWER_ROOT');
        assert(isfile(fullfile(matpowerRoot,'lib','runpf.m')), ...
            'FSA:MATPOWER', 'Set MATPOWER_ROOT to the MATPOWER directory.');
        addpath(genpath(matpowerRoot), '-end');
    end
    algorithmDirectory = fullfile(root,'algorithm','main');
    algorithmName = 'OPFPI_M3E';
    parameters = {};
    if isRWMOP
        problemDirectory = fullfile(platemoRoot,'Problems','Multi-objective optimization','RWMOPs');
        problemName = char(caseName);
    elseif isIEEE118
        algorithmDirectory = fullfile(root,'algorithm','IEEE118');
        algorithmName = 'FSA_CMOEA_IEEE118';
        problemDirectory = fullfile(root,'problems','IEEE118');
        problemName = 'IEEE118_MOOPF_Scalability';
    elseif caseIndex<=4
        problemDirectory = fullfile(root,'problems','IEEE30');
        problemName = 'IEEE30_MOOPF_LEGACY';
        parameters = {'parameter',{caseIndex,0}};
    else
        problemDirectory = fullfile(root,'problems','IEEE57');
        problemName = 'IEEE57_MOOPF_LEGACY';
        parameters = {'parameter',{caseIndex-4,0}};
    end
    addpath(problemDirectory,algorithmDirectory,'-begin');
    assert(strcmp(fileparts(which(algorithmName)),algorithmDirectory) && ...
        strcmp(fileparts(which(problemName)),problemDirectory) && ...
        strcmp(fileparts(which('opfpim3e_problem_profile')),algorithmDirectory), ...
        'FSA:PathConflict', 'A source file is shadowed. Start a fresh MATLAB session.');
    rng(options.Seed,'twister');
    problemFactory = str2func(problemName);
    problem = problemFactory('N',100,'maxFE',options.MaxFE,parameters{:});
    algorithmFactory = str2func(algorithmName);
    algorithm = algorithmFactory('save',0,'run',1,'outputFcn',@(varargin)[]);
    algorithm.Solve(problem);
    population = algorithm.result{end,2};
    result = struct('Case',caseName,'Seed',options.Seed,'MaxFE',options.MaxFE, ...
        'Evaluations',problem.FE,'Dec',population.decs, ...
        'Obj',population.objs,'Con',population.cons);
    if isprop(algorithm,'diagnostics')
        result.Diagnostics = algorithm.diagnostics;
    end
end

