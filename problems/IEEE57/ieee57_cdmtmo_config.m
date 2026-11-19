function cfg = ieee57_cdmtmo_config()
%IEEE57_CDMTMO_CONFIG IEEE 57节点四目标MOOPF工程基准配置（V2）
%
% 统一模型口径：
%   - 基础网络：MATPOWER case57；
%   - 33个控制变量：6 PG + 7 VG + 17 TAP + 3 QC；
%   - PG边界：case57 PMIN/PMAX；
%   - VG边界：发电机母线VMIN/VMAX；
%   - TAP：[0.90,1.10]，步长0.02；
%   - QC：[0,20] MVAr，整数，覆盖原始BS；
%   - 状态约束：7台发电机PG/QG及50个PQ母线电压；
%   - 原case57没有可靠RATE_A，故不构造线路热限约束；
%   - 主要评价指标：HV、IGD。
%
% 本配置以项目自身统一、工程合理、可复现和公平为原则，
% 不要求机械复现某一篇论文的全部参数。

%% 1. 环境检查
requiredNames = {'loadcase','case57','idx_bus','idx_gen','idx_brch'};
for k = 1:numel(requiredNames)
    if isempty(which(requiredNames{k}))
        error('未找到MATPOWER函数或算例：%s',requiredNames{k});
    end
end

[~,~,~,~,BUS_I_COL,BUS_TYPE_COL,~,~,~,BS_COL,~,~,~,~,~,VMAX_COL,VMIN_COL] = idx_bus;
[GEN_BUS_COL,~,~,QMAX_COL,QMIN_COL,VG_COL,~,GEN_STATUS_COL,PMAX_COL,PMIN_COL] = idx_gen;
[~,~,~,~,~,RATE_A_COL,~,~,TAP_COL,SHIFT_COL,BR_STATUS_COL] = idx_brch;

baseMpc = loadcase('case57');

%% 2. 元数据
cfg.meta.id = 'IEEE57_MOOPF_ENGINEERING_V2';
cfg.meta.version = 'CONFIG_ENGINEERING_V2';
cfg.meta.status = 'MODEL_SPEC_V2_LOCKED';
cfg.meta.base_case = 'case57';
cfg.meta.matpower_version_expected = '8.1';
cfg.meta.benchmark_mode = 'engineering_consistent_benchmark';
cfg.meta.exact_paper_reproduction_required = false;
cfg.meta.qc_mode = 'replace_original_BS';
cfg.meta.enforce_q_limits_in_pf = false;
cfg.meta.line_limit_policy = 'OMIT_UNDOCUMENTED_LINE_THERMAL_LIMITS';
cfg.meta.primary_metrics = {'HV','IGD'};

cfg.baseMVA = baseMpc.baseMVA;
cfg.number_of_buses = size(baseMpc.bus,1);
cfg.number_of_generators = size(baseMpc.gen,1);
cfg.number_of_branches = size(baseMpc.branch,1);
cfg.D = 33;
cfg.M_available = 4;
cfg.objective_names = {'FuelCost';'ActivePowerLoss';'Emission';'LoadVoltageDeviation'};

if cfg.number_of_buses ~= 57 || cfg.number_of_generators ~= 7 || cfg.number_of_branches ~= 80
    error('当前case57结构与冻结模型不一致。');
end

%% 3. 控制变量分组
cfg.index.PG  = 1:6;
cfg.index.VG  = 7:13;
cfg.index.TAP = 14:30;
cfg.index.QC  = 31:33;

cfg.mapping.pg_buses = [2;3;6;8;9;12];
cfg.mapping.vg_buses = [1;2;3;6;8;9;12];
cfg.mapping.tap_rows = [19;20;31;35;36;37;41;46;54;58;59;65;66;71;73;76;80];
cfg.mapping.qc_buses = [18;25;53];
cfg.mapping.slack_bus = 1;

onlineGen = baseMpc.gen(:,GEN_STATUS_COL) > 0;
[foundPG,cfg.mapping.pg_gen_rows] = ismember(cfg.mapping.pg_buses,baseMpc.gen(:,GEN_BUS_COL));
[foundVGGen,cfg.mapping.vg_gen_rows] = ismember(cfg.mapping.vg_buses,baseMpc.gen(:,GEN_BUS_COL));
[foundVGBus,cfg.mapping.vg_bus_rows] = ismember(cfg.mapping.vg_buses,baseMpc.bus(:,BUS_I_COL));
[foundQC,cfg.mapping.qc_bus_rows] = ismember(cfg.mapping.qc_buses,baseMpc.bus(:,BUS_I_COL));

if ~all(foundPG) || ~all(onlineGen(cfg.mapping.pg_gen_rows))
    error('PG控制变量到在线发电机的映射失败。');
end
if ~all(foundVGGen) || ~all(onlineGen(cfg.mapping.vg_gen_rows))
    error('VG控制变量到在线发电机的映射失败。');
end
if ~all(foundVGBus)
    error('VG控制变量到母线的映射失败。');
end
if ~all(foundQC)
    error('QC控制变量到母线的映射失败。');
end

cfg.mapping.slack_gen_row = find(baseMpc.gen(:,GEN_BUS_COL)==cfg.mapping.slack_bus & onlineGen,1);
if isempty(cfg.mapping.slack_gen_row)
    error('未找到平衡母线在线发电机。');
end

cfg.mapping.pq_bus_rows = find(baseMpc.bus(:,BUS_TYPE_COL)==1);
cfg.mapping.pq_buses = baseMpc.bus(cfg.mapping.pq_bus_rows,BUS_I_COL);

activeBranch = baseMpc.branch(:,BR_STATUS_COL)>0;
derivedTapRows = find(activeBranch & (abs(baseMpc.branch(:,TAP_COL))>eps | abs(baseMpc.branch(:,SHIFT_COL))>eps));
if ~isequal(derivedTapRows(:),cfg.mapping.tap_rows(:))
    error('冻结的17条变压器支路与当前case57不一致。');
end

%% 4. 变量名称、边界、类型和步长
names = strings(cfg.D,1);
for k = 1:numel(cfg.mapping.pg_buses), names(cfg.index.PG(k)) = "PG"+string(cfg.mapping.pg_buses(k)); end
for k = 1:numel(cfg.mapping.vg_buses), names(cfg.index.VG(k)) = "VG"+string(cfg.mapping.vg_buses(k)); end
for k = 1:numel(cfg.mapping.tap_rows), names(cfg.index.TAP(k)) = "T"+string(cfg.mapping.tap_rows(k)); end
for k = 1:numel(cfg.mapping.qc_buses), names(cfg.index.QC(k)) = "QC"+string(cfg.mapping.qc_buses(k)); end
cfg.variable.names = names;

pgLower = baseMpc.gen(cfg.mapping.pg_gen_rows,PMIN_COL);
pgUpper = baseMpc.gen(cfg.mapping.pg_gen_rows,PMAX_COL);
vgLower = baseMpc.bus(cfg.mapping.vg_bus_rows,VMIN_COL);
vgUpper = baseMpc.bus(cfg.mapping.vg_bus_rows,VMAX_COL);

cfg.variable.lower = [pgLower;vgLower;0.90*ones(17,1);zeros(3,1)];
cfg.variable.upper = [pgUpper;vgUpper;1.10*ones(17,1);20*ones(3,1)];
cfg.variable.kind = [repmat("continuous",13,1);repmat("discrete_step",17,1);repmat("integer",3,1)];
cfg.variable.step = [nan(13,1);repmat(0.02,17,1);ones(3,1)];
cfg.variable.encoding = 'real_with_post_decoding';

%% 5. 状态约束
cfg.constraint.generator_buses = cfg.mapping.vg_buses;
cfg.constraint.generator_rows = cfg.mapping.vg_gen_rows;
cfg.constraint.pg_min = baseMpc.gen(cfg.mapping.vg_gen_rows,PMIN_COL);
cfg.constraint.pg_max = baseMpc.gen(cfg.mapping.vg_gen_rows,PMAX_COL);
cfg.constraint.qg_min = baseMpc.gen(cfg.mapping.vg_gen_rows,QMIN_COL);
cfg.constraint.qg_max = baseMpc.gen(cfg.mapping.vg_gen_rows,QMAX_COL);
cfg.constraint.load_voltage_min = baseMpc.bus(cfg.mapping.pq_bus_rows,VMIN_COL);
cfg.constraint.load_voltage_max = baseMpc.bus(cfg.mapping.pq_bus_rows,VMAX_COL);
cfg.constraint.pg_scale = max(cfg.constraint.pg_max-cfg.constraint.pg_min,1);
cfg.constraint.qg_scale = max(cfg.constraint.qg_max-cfg.constraint.qg_min,1);
cfg.constraint.voltage_scale = max(cfg.constraint.load_voltage_max-cfg.constraint.load_voltage_min,eps);
cfg.constraint.power_flow_failure_penalty = 1e6;
cfg.constraint.feasibility_tolerance = 1e-6;
cfg.constraint.normalized = true;
cfg.constraint.expected_count = 128;
cfg.constraint.include_line_limits = false;

%% 6. 燃料成本系数（精确采用case57.gencost）
if size(baseMpc.gencost,1) ~= cfg.number_of_generators || any(baseMpc.gencost(:,1)~=2) || any(baseMpc.gencost(:,4)~=3)
    error('case57.gencost不是预期的7行二次多项式模型。');
end
cfg.cost.a = baseMpc.gencost(:,7);
cfg.cost.b = baseMpc.gencost(:,6);
cfg.cost.c = baseMpc.gencost(:,5);
cfg.cost.source = 'MATPOWER case57.gencost';
cfg.cost.pg_unit = 'MW';

%% 7. 排放系数
% E_i = alpha_i + beta_i*p_i + gamma_i*p_i^2 + xi_i*exp(lambda_i*p_i)
% p_i = PG_i/baseMVA
cfg.emission.alpha = [0.040;0.030;0.040;0.035;0.050;0.045;0.060];
cfg.emission.beta  = [-0.050;-0.060;-0.050;-0.030;-0.050;-0.040;-0.050];
cfg.emission.gamma = [0.060;0.050;0.040;0.035;0.045;0.050;0.050];
cfg.emission.xi = [0.00002;0.00005;0.00001;0.00002;0.00004;0.00001;0.00001];
cfg.emission.lambda = [0.5;1.5;1.0;0.5;2.0;2.0;1.5];
cfg.emission.pg_unit = 'per_unit_on_baseMVA';
cfg.emission.source = 'fixed_project_benchmark_table';

%% 8. 线路容量策略
cfg.line.rateA = baseMpc.branch(:,RATE_A_COL);
cfg.line.include_constraints = false;
cfg.line.source = 'case57 RATE_A is zero/undocumented; no artificial limits are introduced';
if any(cfg.line.rateA(activeBranch)>0)
    error('当前case57出现非零RATE_A，与V2冻结口径不一致。');
end

%% 9. 目标组合
cfg.case.names = ["Case 1";"Case 2";"Case 3";"Case 4"];
cfg.case.objective_mask = logical([1 1 0 0;1 0 1 0;1 1 1 0;1 1 1 1]);
cfg.case.objective_description = ["FuelCost + ActivePowerLoss";"FuelCost + Emission";"FuelCost + ActivePowerLoss + Emission";"FuelCost + ActivePowerLoss + Emission + LoadVoltageDeviation"];

%% 10. 可选诊断可行点（正式实验默认不使用）
cfg.reference.feasible_seed_available = false;
cfg.reference.feasible_seed_dec = [];
cfg.reference.feasible_seed_file = '';
configFolder = fileparts(mfilename('fullpath'));
seedFile = fullfile(configFolder,'CompatibilityAudit','IEEE57_V2_FEASIBILITY_SEARCH','IEEE57_V2_STRICTLY_FEASIBLE_DIAGNOSTIC_POINT.mat');
if exist(seedFile,'file')
    S = load(seedFile,'diagnosticPoint');
    if isfield(S,'diagnosticPoint') && isfield(S.diagnosticPoint,'decision') && ...
            numel(S.diagnosticPoint.decision)==cfg.D && ...
            isfield(S.diagnosticPoint,'strictFeasible') && logical(S.diagnosticPoint.strictFeasible)
        cfg.reference.feasible_seed_available = true;
        cfg.reference.feasible_seed_dec = reshape(double(S.diagnosticPoint.decision),1,cfg.D);
        cfg.reference.feasible_seed_file = seedFile;
    end
end
cfg.reference.formal_use_allowed = false;

%% 11. 最终一致性检查
if numel(cfg.variable.names)~=cfg.D || numel(cfg.variable.lower)~=cfg.D || ...
        numel(cfg.variable.upper)~=cfg.D || numel(cfg.variable.kind)~=cfg.D || numel(cfg.variable.step)~=cfg.D
    error('变量配置长度与D不一致。');
end
if any(~isfinite(cfg.variable.lower)) || any(~isfinite(cfg.variable.upper)) || any(cfg.variable.lower>cfg.variable.upper)
    error('变量边界无效。');
end
if numel(cfg.mapping.pq_bus_rows)~=50 || numel(cfg.mapping.pg_gen_rows)~=6 || ...
        numel(cfg.mapping.vg_gen_rows)~=7 || numel(cfg.mapping.tap_rows)~=17 || numel(cfg.mapping.qc_bus_rows)~=3
    error('IEEE57控制变量或PQ母线数量与冻结模型不一致。');
end
computedConstraintCount = 4*numel(cfg.constraint.generator_rows)+2*numel(cfg.mapping.pq_bus_rows);
if computedConstraintCount ~= cfg.constraint.expected_count
    error('约束数量应为%d，当前为%d。',cfg.constraint.expected_count,computedConstraintCount);
end
end
