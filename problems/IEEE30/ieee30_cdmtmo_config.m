
function cfg = ieee30_cdmtmo_config()
%IEEE30_CDMTMO_CONFIG IEEE 30节点多目标最优潮流统一配置
%
% 本配置用于建立数据来源明确、计算稳定、约束自洽的
% IEEE 30节点多目标最优潮流工程基准。
%
% 返回结构体cfg，统一保存：
%   1. 24个控制变量的顺序、类型、步长和上下界；
%   2. 控制变量与MATPOWER节点、发电机和支路的映射；
%   3. 发电机有功、无功和负荷节点电压约束；
%   4. 燃料成本和污染排放系数；
%   5. 工程线路容量；
%   6. 论文Table 3的四组舍入参考解；
%   7. 已完成的模型审计信息。
%
% 重要建模决定：
%   - 基础网络使用MATPOWER case_ieee30；
%   - 线路容量使用MATPOWER case30中的RATE_A，并按支路顺序映射；
%   - QC直接覆盖case_ieee30原始BS，不进行叠加；
%   - 潮流计算不启用MATPOWER无功限自动修复；
%   - 变压器分接比范围为0.90到1.10，步长暂按0.02处理；
%   - 并联无功补偿变量为0到5 MVAr之间的整数；
%   - 论文结果仅作为合理性参照，不要求逐位精确复现。
%
% 本函数只读取MATPOWER数据，不修改任何MATPOWER原始文件。

%% 1. 环境检查

requiredNames = {
    'loadcase'
    'define_constants'
    'case_ieee30'
    'case30'
    };

for k = 1:numel(requiredNames)
    if isempty(which(requiredNames{k}))
        error('未找到MATPOWER函数或算例：%s', requiredNames{k});
    end
end

define_constants;

baseMpc = loadcase('case_ieee30');
rateMpc = loadcase('case30');

%% 2. 基本元数据

cfg.meta.id = 'IEEE30_MOOPF_ENGINEERING_V1';
cfg.meta.version = 'CONFIG_ENGINEERING_V1';

cfg.meta.description = ...
    ['IEEE 30-bus multi-objective optimal power flow benchmark ' ...
     'based on MATPOWER case_ieee30, with explicitly defined ' ...
     'objectives, controls, constraints and engineering line limits.'];

cfg.meta.created_from = ...
    ['MATPOWER 8.1 case_ieee30, MATPOWER case30 RATE_A, ' ...
     'and CDMTMO paper reference parameters'];

cfg.meta.base_case = 'case_ieee30';
cfg.meta.line_limit_case = 'case30';
cfg.meta.matpower_version_expected = '8.1';

cfg.meta.benchmark_mode = ...
    'engineering_consistent_benchmark';

cfg.meta.exact_paper_reproduction_required = false;
cfg.meta.paper_exact_source_recovered = false;
cfg.meta.table3_is_rounded_reference = true;

cfg.meta.qc_mode = 'replace_original_BS';
cfg.meta.enforce_q_limits_in_pf = false;

cfg.baseMVA = baseMpc.baseMVA;

cfg.number_of_buses = size(baseMpc.bus,1);
cfg.number_of_generators = size(baseMpc.gen,1);
cfg.number_of_branches = size(baseMpc.branch,1);

cfg.D = 24;
cfg.M_available = 4;

cfg.objective_names = {
    'FuelCost'
    'ActivePowerLoss'
    'Emission'
    'LoadVoltageDeviation'
    };

%% 3. 控制变量分组

cfg.index.PG  = 1:5;
cfg.index.VG  = 6:11;
cfg.index.TAP = 12:15;
cfg.index.QC  = 16:24;

% 变量顺序：
% PG2, PG5, PG8, PG11, PG13,
% VG1, VG2, VG5, VG8, VG11, VG13,
% T11, T12, T15, T36,
% QC10, QC12, QC15, QC17, QC20, QC21, QC23, QC24, QC29

cfg.variable.names = [
    "PG2"
    "PG5"
    "PG8"
    "PG11"
    "PG13"
    "VG1"
    "VG2"
    "VG5"
    "VG8"
    "VG11"
    "VG13"
    "T11"
    "T12"
    "T15"
    "T36"
    "QC10"
    "QC12"
    "QC15"
    "QC17"
    "QC20"
    "QC21"
    "QC23"
    "QC24"
    "QC29"
    ];

% 注意：lower和upper均保持为24×1列向量

cfg.variable.lower = [
    20
    15
    10
    10
    12
    0.95
    0.95
    0.95
    0.95
    0.95
    0.95
    0.90
    0.90
    0.90
    0.90
    0
    0
    0
    0
    0
    0
    0
    0
    0
    ];

cfg.variable.upper = [
    80
    50
    35
    30
    40
    1.10
    1.10
    1.10
    1.10
    1.10
    1.10
    1.10
    1.10
    1.10
    1.10
    5
    5
    5
    5
    5
    5
    5
    5
    5
    ];

cfg.variable.kind = [
    repmat("continuous",1,11), ...
    repmat("discrete_step",1,4), ...
    repmat("integer",1,9)
    ]';

cfg.variable.step = [
    nan(1,11), ...
    repmat(0.02,1,4), ...
    ones(1,9)
    ]';

cfg.variable.encoding = 'real_with_post_decoding';

%% 4. 控制变量到MATPOWER的映射

cfg.mapping.pg_buses = [2; 5; 8; 11; 13];
cfg.mapping.vg_buses = [1; 2; 5; 8; 11; 13];
cfg.mapping.tap_rows = [11; 12; 15; 36];
cfg.mapping.qc_buses = [10; 12; 15; 17; 20; 21; 23; 24; 29];

[foundPG,cfg.mapping.pg_gen_rows] = ismember( ...
    cfg.mapping.pg_buses, ...
    baseMpc.gen(:,GEN_BUS));

[foundVGGen,cfg.mapping.vg_gen_rows] = ismember( ...
    cfg.mapping.vg_buses, ...
    baseMpc.gen(:,GEN_BUS));

[foundVGBus,cfg.mapping.vg_bus_rows] = ismember( ...
    cfg.mapping.vg_buses, ...
    baseMpc.bus(:,BUS_I));

[foundQC,cfg.mapping.qc_bus_rows] = ismember( ...
    cfg.mapping.qc_buses, ...
    baseMpc.bus(:,BUS_I));

if ~all(foundPG)
    error('PG节点到发电机表的映射失败。');
end

if ~all(foundVGGen)
    error('VG节点到发电机表的映射失败。');
end

if ~all(foundVGBus)
    error('VG节点到节点表的映射失败。');
end

if ~all(foundQC)
    error('QC节点到节点表的映射失败。');
end

if any(cfg.mapping.tap_rows < 1) || ...
        any(cfg.mapping.tap_rows > size(baseMpc.branch,1))
    error('可调变压器支路序号超出基础算例范围。');
end

cfg.mapping.slack_bus = 1;

cfg.mapping.slack_gen_row = find( ...
    baseMpc.gen(:,GEN_BUS) == cfg.mapping.slack_bus, ...
    1);

if isempty(cfg.mapping.slack_gen_row)
    error('未找到平衡节点发电机。');
end

cfg.mapping.pq_bus_rows = find( ...
    baseMpc.bus(:,BUS_TYPE) == PQ);

cfg.mapping.pq_buses = ...
    baseMpc.bus(cfg.mapping.pq_bus_rows,BUS_I);

%% 5. 发电机状态约束

% 发电机顺序：
% 节点1、2、5、8、11、13

cfg.constraint.generator_buses = [
    1
    2
    5
    8
    11
    13
    ];

cfg.constraint.pg_min = [
    50
    20
    15
    10
    10
    12
    ];

cfg.constraint.pg_max = [
    250
    80
    50
    35
    30
    40
    ];

cfg.constraint.qg_min = [
    -20
    -20
    -15
    -15
    -10
    -15
    ];

cfg.constraint.qg_max = [
    150
    60
    62.5
    48.7
    40
    44.7
    ];

cfg.constraint.load_voltage_min = 0.95;
cfg.constraint.load_voltage_max = 1.05;

% 潮流失败时使用的大惩罚值
cfg.constraint.power_flow_failure_penalty = 1e6;

% 数值可行性判定容差
cfg.constraint.feasibility_tolerance = 1e-8;

%% 6. 发电成本系数

% Cost_i = a_i + b_i*PG_i + c_i*PG_i^2

cfg.cost.a = zeros(6,1);

cfg.cost.b = [
    2.00
    1.75
    1.00
    3.25
    3.00
    3.00
    ];

cfg.cost.c = [
    0.00375
    0.01750
    0.06250
    0.00834
    0.02500
    0.02500
    ];

cfg.cost.pg_unit = 'MW';
cfg.cost.output_unit = '$/h';

%% 7. 排放系数

% 计算口径：
%
% p_i = PG_i / baseMVA
%
% Emission_i =
%   0.01*(alpha_i + beta_i*p_i + gamma_i*p_i^2)
%   + xi_i*exp(lambda_i*p_i)

cfg.emission.alpha = [
    4.091
    2.543
    4.258
    5.326
    4.258
    6.131
    ];

cfg.emission.beta = [
   -5.554
   -6.074
   -5.094
   -3.550
   -5.094
   -5.555
    ];

cfg.emission.gamma = [
    6.490
    5.638
    4.586
    3.380
    4.586
    5.151
    ];

cfg.emission.xi = [
    0.000200
    0.000500
    0.000001
    0.002000
    0.000001
    0.000010
    ];

cfg.emission.lambda = [
    2.857
    3.333
    8.000
    2.000
    8.000
    6.667
    ];

cfg.emission.pg_unit = 'per_unit_on_baseMVA';
cfg.emission.output_unit = 'ton/h';

%% 8. 工程线路容量

% case_ieee30的RATE_A全部为0。
% 本基准采用case30中的RATE_A，并按照完全一致的支路起止节点
% 和支路顺序映射到case_ieee30。

sameTopology = ...
    all(baseMpc.branch(:,F_BUS) == ...
        rateMpc.branch(:,F_BUS)) && ...
    all(baseMpc.branch(:,T_BUS) == ...
        rateMpc.branch(:,T_BUS));

if ~sameTopology
    error(['case30与case_ieee30的支路起止节点或' ...
           '支路顺序不一致，不能映射RATE_A。']);
end

cfg.line.rateA = rateMpc.branch(:,RATE_A);

if any(~isfinite(cfg.line.rateA))
    error('工程RATE_A中存在非有限数值。');
end

if any(cfg.line.rateA <= 0)
    error('工程RATE_A中存在非正值。');
end

% 保留旧字段名，兼容前期审计脚本
cfg.line.rateA_candidate = cfg.line.rateA;

cfg.line.source = ...
    ['MATPOWER case30 RATE_A mapped to case_ieee30 ' ...
     'by identical branch endpoints and branch order'];

cfg.line.source_traceable = true;
cfg.line.topology_mapping_verified = true;

% 是否确认与CDMTMO论文原始Smax完全相同
cfg.line.paper_source_verified = false;

% 当前模型不定位为严格论文复现
cfg.line.use_in_formal_reproduction = false;

% 在本工程一致性基准中正式启用线路容量
cfg.line.use_in_formal_benchmark = true;

% 保留前期兼容字段
cfg.line.use_in_candidate_engineering_test = true;

cfg.line.audit_max_loading_percent_cases_1_to_4 = [
    58.949349
    65.054780
    59.738551
    72.532242
    ];

cfg.line.audit_violating_count_cases_1_to_4 = [
    0
    0
    0
    0
    ];

%% 9. 四个参考目标组合

% 可用目标顺序：
% 1 FuelCost
% 2 ActivePowerLoss
% 3 Emission
% 4 LoadVoltageDeviation

cfg.case.names = [
    "Case 1"
    "Case 2"
    "Case 3"
    "Case 4"
    ];

cfg.case.objective_mask = logical([
    1 1 0 0
    1 0 1 0
    1 1 1 0
    1 1 1 1
    ]);

cfg.case.objective_description = [
    "FuelCost + ActivePowerLoss"
    "FuelCost + Emission"
    "FuelCost + ActivePowerLoss + Emission"
    "FuelCost + ActivePowerLoss + Emission + LoadVoltageDeviation"
    ];

%% 10. 论文Table 3中的舍入参考控制变量

% 以下数据只作为模型合理性校验和回归测试参考，
% 不作为必须逐位复现的正式优化结果。

pgReference = [
    52.4349    58.2646    59.9079    60.9826
    34.0744    27.1172    38.4653    37.0260
    35.0000    34.9958    34.9983    34.9663
    29.8462    26.1937    29.9225    29.9981
    21.9436    25.2784    32.1052    30.1204
    ];

vgReference = [
    1.0718    1.0704    1.0662    1.0282
    1.0600    1.0584    1.0574    1.0211
    1.0352    1.0330    1.0350    1.0004
    1.0425    1.0397    1.0424    1.0059
    1.0847    1.0768    1.0707    1.0050
    1.0480    1.0611    1.0577    1.0200
    ];

tapReference = [
    1.06    1.06    1.04    1.02
    0.94    0.92    0.94    0.90
    0.98    1.00    1.00    1.00
    0.98    0.98    0.98    0.96
    ];

qcReference = [
    5    3    4    4
    3    3    3    2
    4    2    3    3
    5    5    5    1
    4    4    4    5
    5    5    5    5
    3    4    4    5
    5    5    5    5
    3    2    3    2
    ];

% 每一行对应一个Case，每一列对应一个控制变量
% 维度：4 × 24
cfg.reference.table3_x = ...
    [pgReference' vgReference' tapReference' qcReference'];

cfg.reference.paper_pg1 = [
    115.2076
    117.1291
     92.2559
     95.1246
    ];

% 列顺序：
% Cost, Ploss, Emission, LVD

cfg.reference.paper_objectives = [
    840.3490    5.1067    0.2432    0.9110
    831.6472    5.5788    0.2476    0.8541
    873.2087    4.2551    0.2199    0.8925
    868.0614    4.8180    0.2230    0.1161
    ];

cfg.reference.is_exact_solution = false;
cfg.reference.purpose = ...
    'rounded_reference_for_mapping_and_regression_audit';

%% 11. 已完成的审计结果

cfg.audit.power_flow_success_cases_1_to_4 = ...
    true(4,1);

cfg.audit.calculated_pg1_cases_1_to_4 = [
    115.179812
    117.100693
     92.231635
     95.083106
    ];

cfg.audit.calculated_objectives_cases_1_to_4 = [
    840.269262    5.078912    0.243045    0.935193
    831.565438    5.550393    0.247364    0.879052
    873.143450    4.230835    0.219761    0.916834
    867.948822    4.776506    0.222766    0.119675
    ];

cfg.audit.pg_violation_cases_1_to_4 = ...
    zeros(4,1);

cfg.audit.qg_violation_cases_1_to_4 = ...
    zeros(4,1);

cfg.audit.voltage_violation_cases_1_to_4 = [
    0.001879241218
    0.001827143883
    0.001601945845
    0
    ];

cfg.audit.line_violation_cases_1_to_4 = ...
    zeros(4,1);

cfg.audit.strict_feasible_from_rounded_table = [
    false
    false
    false
    true
    ];

cfg.audit.note = [
    "All four rounded reference controls produce converged power-flow solutions."
    "PG and QG constraints are satisfied for all four rounded reference controls."
    "Cases 1-3 show only small voltage exceedances caused by rounded controls."
    "Case 4 is strictly feasible using the rounded reference controls."
    "Engineering RATE_A limits are inactive at all four reference points."
    "The benchmark does not require exact numerical reproduction of the paper."
    ];

%% 12. 最终一致性检查

if numel(cfg.variable.names) ~= cfg.D
    error('变量名称数量与D不一致。');
end

if numel(cfg.variable.lower) ~= cfg.D
    error('变量下界数量与D不一致。');
end

if numel(cfg.variable.upper) ~= cfg.D
    error('变量上界数量与D不一致。');
end

if numel(cfg.variable.kind) ~= cfg.D
    error('变量类型数量与D不一致。');
end

if numel(cfg.variable.step) ~= cfg.D
    error('变量步长数量与D不一致。');
end

if ~isequal(size(cfg.variable.names),[cfg.D,1])
    error('变量名称必须为24×1列向量。');
end

if ~isequal(size(cfg.variable.lower),[cfg.D,1])
    error('变量下界必须为24×1列向量。');
end

if ~isequal(size(cfg.variable.upper),[cfg.D,1])
    error('变量上界必须为24×1列向量。');
end

if ~isequal(size(cfg.variable.kind),[cfg.D,1])
    error('变量类型必须为24×1列向量。');
end

if ~isequal(size(cfg.variable.step),[cfg.D,1])
    error('变量步长必须为24×1列向量。');
end

if any(~isfinite(cfg.variable.lower)) || ...
        any(~isfinite(cfg.variable.upper))
    error('变量边界中存在非有限数值。');
end

if any(cfg.variable.lower > cfg.variable.upper)
    error('存在变量下界大于上界的情况。');
end

if ~isequal(size(cfg.reference.table3_x),[4,cfg.D])
    error('Table 3参考控制变量矩阵必须为4×24。');
end

if size(cfg.reference.paper_objectives,1) ~= 4 || ...
        size(cfg.reference.paper_objectives,2) ~= cfg.M_available
    error('论文参考目标矩阵尺寸应为4×4。');
end

if numel(cfg.mapping.pq_buses) ~= 24
    error('负荷节点数量不是24。');
end

if numel(cfg.mapping.pg_gen_rows) ~= 5
    error('非平衡发电机映射数量不是5。');
end

if numel(cfg.mapping.vg_gen_rows) ~= 6
    error('发电机电压映射数量不是6。');
end

if numel(cfg.mapping.tap_rows) ~= 4
    error('可调变压器数量不是4。');
end

if numel(cfg.mapping.qc_bus_rows) ~= 9
    error('并联无功补偿节点数量不是9。');
end

if numel(cfg.line.rateA) ~= cfg.number_of_branches
    error('线路容量数量与支路数量不一致。');
end

end

