
function [F,G,detail,Xdecoded] = ieee30_moopf_evaluate_legacy(X,cfg)
%IEEE30_MOOPF_EVALUATE IEEE 30???MOOPF?????????
%
% ?????
%   [F,G] = ieee30_moopf_evaluate(X);
%   [F,G,detail,Xdecoded] = ieee30_moopf_evaluate(X,cfg);
%
% ?????
%   X   - N?24???????????????24???????????
%   cfg - ieee30_cdmtmo_config????????????
%
% ?????
%   F - N?4????????
%       1 FuelCost
%       2 ActivePowerLoss
%       3 Emission
%       4 LoadVoltageDeviation
%
%   G - N?113???????????latEMO????????
%       G <= 0????????
%       G > 0 ????????
%
%       ????????
%       1-6     PG???
%       7-12    PG???
%       13-18   QG???
%       19-24   QG???
%       25-48   ????????????
%       49-72   ????????????
%       73-113  ?????????
%
%   detail   - ?????????
%   Xdecoded - ?????????N?24??????
%
% ????????
%   - ????????ATPOWER case_ieee30
%   - ???????????ATPOWER case30 RATE_A????????
%   - QC??????BS
%   - TAP???0.02
%   - QC?????
%   - ?????ATPOWER???????????
%
% ?????
%   ???????????dx_bus??dx_gen??dx_brch???????????
%   ???_COL????????ATPOWER??D?????latEMO??D.m?????
%
% ??????????????ATPOWER????????

%% 1. ???????

if nargin < 2 || isempty(cfg)
    cfg = ieee30_cdmtmo_config();
end

if ~isstruct(cfg) || ...
        ~isfield(cfg,'D') || ...
        ~isfield(cfg,'mapping') || ...
        ~isfield(cfg,'constraint')
    error('IEEE legacy evaluation compatibility error.');
end

if ~isfield(cfg,'line') || ...
        ~isfield(cfg.line,'rateA') || ...
        ~cfg.line.use_in_formal_benchmark
    error('IEEE legacy evaluation compatibility error.');
end

%% 2. ????????????

[XdecodedRaw,decodeInfo] = ...
    ieee30_moopf_decode(X,cfg);

if isvector(XdecodedRaw)
    Xdecoded = reshape(XdecodedRaw,1,cfg.D);
else
    Xdecoded = XdecodedRaw;
end

numberOfSolutions = size(Xdecoded,1);

%% 3. MATPOWER???????

requiredNames = {
    'loadcase'
    'runpf'
    'mpoption'
    'idx_bus'
    'idx_gen'
    'idx_brch'
    'case_ieee30'
    };

for k = 1:numel(requiredNames)
    if isempty(which(requiredNames{k}))
        error('IEEE legacy evaluation compatibility error.');
    end
end

%% 4. ??????MATPOWER???
%
% ?????efine_constants???PD?????????
% ?????latEMO???PD.m??????????????

[ ...
    PQ_TYPE, ...
    PV_TYPE, ...
    REF_TYPE, ...
    NONE_TYPE, ...
    BUS_I_COL, ...
    BUS_TYPE_COL, ...
    PD_COL, ...
    QD_COL, ...
    GS_COL, ...
    BS_COL, ...
    BUS_AREA_COL, ...
    VM_COL, ...
    VA_COL, ...
    BASE_KV_COL, ...
    ZONE_COL, ...
    VMAX_COL, ...
    VMIN_COL ...
    ] = idx_bus;

[ ...
    GEN_BUS_COL, ...
    PG_COL, ...
    QG_COL, ...
    QMAX_COL, ...
    QMIN_COL, ...
    VG_COL, ...
    MBASE_COL, ...
    GEN_STATUS_COL, ...
    PMAX_COL, ...
    PMIN_COL ...
    ] = idx_gen;

[ ...
    F_BUS_COL, ...
    T_BUS_COL, ...
    BR_R_COL, ...
    BR_X_COL, ...
    BR_B_COL, ...
    RATE_A_COL, ...
    RATE_B_COL, ...
    RATE_C_COL, ...
    TAP_COL, ...
    SHIFT_COL, ...
    BR_STATUS_COL, ...
    PF_COL, ...
    QF_COL, ...
    PT_COL, ...
    QT_COL ...
    ] = idx_brch;

% ???????????????????????????????????
% ?????COL?????????????????????????????

%% 5. ?????????

baseMpc = loadcase(cfg.meta.base_case);

if size(baseMpc.bus,1) ~= cfg.number_of_buses
    error('IEEE legacy evaluation compatibility error.');
end

if size(baseMpc.gen,1) ~= cfg.number_of_generators
    error('IEEE legacy evaluation compatibility error.');
end

if size(baseMpc.branch,1) ~= cfg.number_of_branches
    error('IEEE legacy evaluation compatibility error.');
end

%% 6. ???????????

[foundStateGen,stateGenRows] = ismember( ...
    cfg.constraint.generator_buses, ...
    baseMpc.gen(:,GEN_BUS_COL));

if ~all(foundStateGen)
    error('IEEE legacy evaluation compatibility error.');
end

numberOfGenerators = numel(stateGenRows);
numberOfPQBuses = numel(cfg.mapping.pq_bus_rows);
numberOfBranches = size(baseMpc.branch,1);

numberOfConstraints = ...
    4*numberOfGenerators + ...
    2*numberOfPQBuses + ...
    numberOfBranches;

if numberOfConstraints ~= 113
    error('IEEE legacy evaluation compatibility error.');
end

%% 7. ????????

F = nan(numberOfSolutions,4);
G = nan(numberOfSolutions,numberOfConstraints);

detail.success = false(numberOfSolutions,1);

detail.pg1 = nan(numberOfSolutions,1);
detail.total_pg = nan(numberOfSolutions,1);
detail.total_pd = nan(numberOfSolutions,1);

detail.pg_violation_sum = ...
    nan(numberOfSolutions,1);

detail.qg_violation_sum = ...
    nan(numberOfSolutions,1);

detail.voltage_violation_sum = ...
    nan(numberOfSolutions,1);

detail.line_violation_sum = ...
    nan(numberOfSolutions,1);

detail.max_line_loading_percent = ...
    nan(numberOfSolutions,1);

detail.max_line_loading_branch = ...
    nan(numberOfSolutions,1);

detail.power_flow_error = ...
    strings(numberOfSolutions,1);

%% 8. ??????

mpopt = mpoption( ...
    'verbose',0, ...
    'out.all',0, ...
    'pf.enforce_q_lims',0, ...
    'exp.use_legacy_core',1);

failurePenalty = ...
    cfg.constraint.power_flow_failure_penalty;

rateA = cfg.line.rateA(:);

if numel(rateA) ~= numberOfBranches
    error('IEEE legacy evaluation compatibility error.');
end

if any(~isfinite(rateA)) || any(rateA <= 0)
    error('IEEE legacy evaluation compatibility error.');
end

%% 9. ???????????

for solutionIndex = 1:numberOfSolutions

    x = Xdecoded(solutionIndex,:);

    % ????????????????????????????????????
    mpc = baseMpc;

    %% 9.1 ??????????????????

    mpc.gen(cfg.mapping.pg_gen_rows,PG_COL) = ...
        x(cfg.index.PG)';

    %% 9.2 ????????????????

    mpc.gen(cfg.mapping.vg_gen_rows,VG_COL) = ...
        x(cfg.index.VG)';

    % ???????????????
    mpc.bus(cfg.mapping.vg_bus_rows,VM_COL) = ...
        x(cfg.index.VG)';

    %% 9.3 ????????????

    mpc.branch(cfg.mapping.tap_rows,TAP_COL) = ...
        x(cfg.index.TAP)';

    %% 9.4 ????????????

    % ?????????QC??????BS?????????
    mpc.bus(cfg.mapping.qc_bus_rows,BS_COL) = ...
        x(cfg.index.QC)';

    %% 9.5 ??????

    try
        [result,success] = runpf(mpc,mpopt);
    catch ME
        result = [];
        success = false;

        detail.power_flow_error(solutionIndex) = string(ME.message);

    end

    detail.success(solutionIndex) = logical(success);

    %% 9.6 ?????????

    if ~success
        F(solutionIndex,:) = failurePenalty;
        G(solutionIndex,:) = failurePenalty;
        continue;
    end

    %% 9.7 ??????????

    pg = result.gen(stateGenRows,PG_COL);
    qg = result.gen(stateGenRows,QG_COL);

    pqVoltage = ...
        result.bus(cfg.mapping.pq_bus_rows,VM_COL);

    %% 9.8 ???1????????

    fuelCost = sum( ...
        cfg.cost.a + ...
        cfg.cost.b .* pg + ...
        cfg.cost.c .* pg.^2);

    %% 9.9 ???2???????????

    onlineGeneratorRows = ...
        result.gen(:,GEN_STATUS_COL) > 0;

    totalPG = sum( ...
        result.gen(onlineGeneratorRows,PG_COL));

    totalPD = sum( ...
        result.bus(:,PD_COL));

    activePowerLoss = ...
        totalPG-totalPD;

    %% 9.10 ???3????????

    pgPerUnit = ...
        pg ./ result.baseMVA;

    emission = sum( ...
        0.01 .* ( ...
            cfg.emission.alpha + ...
            cfg.emission.beta .* pgPerUnit + ...
            cfg.emission.gamma .* pgPerUnit.^2) + ...
        cfg.emission.xi .* ...
            exp(cfg.emission.lambda .* pgPerUnit));

    %% 9.11 ???4??????????????

    loadVoltageDeviation = ...
        sum(abs(pqVoltage-1.0));

    F(solutionIndex,:) = [
        fuelCost
        activePowerLoss
        emission
        loadVoltageDeviation
        ]';

    %% 9.12 PG???
    %
    % g <= 0?????????

    gPgLower = ...
        cfg.constraint.pg_min-pg;

    gPgUpper = ...
        pg-cfg.constraint.pg_max;

    %% 9.13 QG???

    gQgLower = ...
        cfg.constraint.qg_min-qg;

    gQgUpper = ...
        qg-cfg.constraint.qg_max;

    %% 9.14 ????????????

    gVoltageLower = ...
        cfg.constraint.load_voltage_min-pqVoltage;

    gVoltageUpper = ...
        pqVoltage-cfg.constraint.load_voltage_max;

    %% 9.15 ????????????

    apparentPowerFrom = hypot( ...
        result.branch(:,PF_COL), ...
        result.branch(:,QF_COL));

    apparentPowerTo = hypot( ...
        result.branch(:,PT_COL), ...
        result.branch(:,QT_COL));

    maximumApparentPower = max( ...
        apparentPowerFrom, ...
        apparentPowerTo);

    gLine = ...
        maximumApparentPower-rateA;

    %% 9.16 ?????????

    G(solutionIndex,:) = [
        gPgLower
        gPgUpper
        gQgLower
        gQgUpper
        gVoltageLower
        gVoltageUpper
        gLine
        ]';

    %% 9.17 ??????

    detail.pg1(solutionIndex) = ...
        result.gen( ...
            cfg.mapping.slack_gen_row, ...
            PG_COL);

    detail.total_pg(solutionIndex) = totalPG;
    detail.total_pd(solutionIndex) = totalPD;

    detail.pg_violation_sum(solutionIndex) = sum( ...
        max(gPgLower,0) + ...
        max(gPgUpper,0));

    detail.qg_violation_sum(solutionIndex) = sum( ...
        max(gQgLower,0) + ...
        max(gQgUpper,0));

    detail.voltage_violation_sum(solutionIndex) = sum( ...
        max(gVoltageLower,0) + ...
        max(gVoltageUpper,0));

    detail.line_violation_sum(solutionIndex) = ...
        sum(max(gLine,0));

    loadingPercent = ...
        100 .* maximumApparentPower ./ rateA;

    [maximumLoading,maximumLoadingBranch] = ...
        max(loadingPercent);

    detail.max_line_loading_percent(solutionIndex) = ...
        maximumLoading;

    detail.max_line_loading_branch(solutionIndex) = ...
        maximumLoadingBranch;
end

%% 10. ?????????

positiveG = max(G,0);

detail.cv_sum = ...
    sum(positiveG,2);

detail.cv_max = ...
    max(positiveG,[],2);

detail.is_feasible = ...
    detail.success & ...
    detail.cv_max <= ...
        cfg.constraint.feasibility_tolerance;

detail.number_of_solutions = ...
    numberOfSolutions;

detail.number_of_objectives = ...
    size(F,2);

detail.number_of_constraints = ...
    size(G,2);

detail.objective_names = ...
    cfg.objective_names;

detail.decode_info = ...
    decodeInfo;

%% 11. ?????????

cursor = 0;

detail.constraint_index.pg_lower = ...
    cursor + (1:numberOfGenerators);

cursor = cursor + numberOfGenerators;

detail.constraint_index.pg_upper = ...
    cursor + (1:numberOfGenerators);

cursor = cursor + numberOfGenerators;

detail.constraint_index.qg_lower = ...
    cursor + (1:numberOfGenerators);

cursor = cursor + numberOfGenerators;

detail.constraint_index.qg_upper = ...
    cursor + (1:numberOfGenerators);

cursor = cursor + numberOfGenerators;

detail.constraint_index.voltage_lower = ...
    cursor + (1:numberOfPQBuses);

cursor = cursor + numberOfPQBuses;

detail.constraint_index.voltage_upper = ...
    cursor + (1:numberOfPQBuses);

cursor = cursor + numberOfPQBuses;

detail.constraint_index.line_upper = ...
    cursor + (1:numberOfBranches);

cursor = cursor + numberOfBranches;

if cursor ~= numberOfConstraints
    error('IEEE legacy evaluation compatibility error.');
end

%% 12. ???????????

detail.index_source = ...
    'MATPOWER idx_bus / idx_gen / idx_brch';

detail.name_collision_avoided = true;

detail.constraint_sign_rule = ...
    'G <= 0 feasible; G > 0 violated';

end


