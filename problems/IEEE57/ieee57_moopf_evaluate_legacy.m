function [F,G,detail,Xdecoded] = ieee57_moopf_evaluate_legacy(X,cfg)
%IEEE57_MOOPF_EVALUATE IEEE57 V2???????28??????????????%
% F???FuelCost, ActivePowerLoss, Emission, LoadVoltageDeviation??% G???????????e-6?????<=0??????
%   1:7    PG?????:14 PG?????%   15:21  QG?????2:28 QG?????%   29:78  PQ????????9:128 PQ????????
if nargin<2 || isempty(cfg), cfg = ieee57_cdmtmo_config(); end
if ~isstruct(cfg) || ~isfield(cfg,'D') || ~isfield(cfg,'mapping') || ~isfield(cfg,'constraint')
    error('IEEE legacy evaluation compatibility error.');
end
[Xraw,decodeInfo] = ieee57_moopf_decode(X,cfg);
if isvector(Xraw), Xdecoded = reshape(Xraw,1,cfg.D); else, Xdecoded = Xraw; end
N = size(Xdecoded,1);

requiredNames = {'loadcase','runpf','mpoption','totcost','idx_bus','idx_gen','idx_brch','case57'};
for k=1:numel(requiredNames)
    if isempty(which(requiredNames{k})), error('IEEE legacy evaluation compatibility error.'); end
end

[~,~,~,~,BUS_I_COL,~,PD_COL,~,~,BS_COL,~,VM_COL,~,~,~,~,~] = idx_bus;
[GEN_BUS_COL,PG_COL,QG_COL,~,~,VG_COL,~,GEN_STATUS_COL,~,~] = idx_gen;
[~,~,~,~,~,~,~,~,TAP_COL,~,BR_STATUS_COL,PF_COL,~,PT_COL,~] = idx_brch;

% MATPOWER 8.1????????????VM=8??F=14??T=16??% ???????????dx_*???????????????????A???????????M??if VM_COL~=8 || PF_COL~=14 || PT_COL~=16
% MATPOWER index sanity check.
if VM_COL~=8 || PF_COL~=14 || PT_COL~=16
    error('IEEE legacy evaluation compatibility error.');
end

baseMpc = loadcase(cfg.meta.base_case);
if size(baseMpc.bus,1)~=cfg.number_of_buses || size(baseMpc.gen,1)~=cfg.number_of_generators || size(baseMpc.branch,1)~=cfg.number_of_branches
    error('IEEE legacy evaluation compatibility error.');
end

stateGenRows = cfg.constraint.generator_rows(:);
numberOfGenerators = numel(stateGenRows);
numberOfPQ = numel(cfg.mapping.pq_bus_rows);
numberOfConstraints = 4*numberOfGenerators+2*numberOfPQ;
if numberOfConstraints~=cfg.constraint.expected_count, error('IEEE legacy evaluation compatibility error.'); end

F = nan(N,4);
G = nan(N,numberOfConstraints);
detail.success = false(N,1);
detail.power_flow_error = strings(N,1);
detail.slack_pg = nan(N,1);
detail.total_pg = nan(N,1);
detail.total_pd = nan(N,1);
detail.minimum_pq_voltage = nan(N,1);
detail.maximum_pq_voltage = nan(N,1);
detail.cost_from_totcost = nan(N,1);
detail.cost_from_coefficients = nan(N,1);
detail.cost_difference = nan(N,1);

mpopt = mpoption('verbose',0,'out.all',0,'pf.alg','NR','pf.tol',1e-8,'pf.nr.max_it',50,'pf.enforce_q_lims',0,'exp.use_legacy_core',1);
penalty = cfg.constraint.power_flow_failure_penalty;
tolerance = cfg.constraint.feasibility_tolerance;
activeBranchRows = find(baseMpc.branch(:,BR_STATUS_COL)>0);

for i=1:N
    x = Xdecoded(i,:);
    mpc = baseMpc;
    mpc.gen(cfg.mapping.pg_gen_rows,PG_COL) = x(cfg.index.PG)';
    mpc.gen(cfg.mapping.vg_gen_rows,VG_COL) = x(cfg.index.VG)';
    mpc.bus(cfg.mapping.vg_bus_rows,VM_COL) = x(cfg.index.VG)';
    mpc.branch(cfg.mapping.tap_rows,TAP_COL) = x(cfg.index.TAP)';
    mpc.bus(cfg.mapping.qc_bus_rows,BS_COL) = x(cfg.index.QC)';
    try
        [result,success] = runpf(mpc,mpopt);
    catch ME
        result=[]; success=false; detail.power_flow_error(i)=string(ME.message);
    end
    detail.success(i)=logical(success);
    if ~success
        F(i,:)=penalty; G(i,:)=penalty; continue;
    end

    pg = result.gen(stateGenRows,PG_COL);
    qg = result.gen(stateGenRows,QG_COL);
    pqVoltage = result.bus(cfg.mapping.pq_bus_rows,VM_COL);

    fuelCost = sum(totcost(result.gencost,result.gen(:,PG_COL)));
    costByCoeff = sum(cfg.cost.a+cfg.cost.b.*pg+cfg.cost.c.*pg.^2);
    activeLoss = sum(result.branch(activeBranchRows,PF_COL)+result.branch(activeBranchRows,PT_COL));
    if activeLoss<0 && activeLoss>-1e-8, activeLoss=0; end
    p = pg./result.baseMVA;
    emission = sum(cfg.emission.alpha+cfg.emission.beta.*p+cfg.emission.gamma.*p.^2+cfg.emission.xi.*exp(cfg.emission.lambda.*p));
    lvd = sum(abs(pqVoltage-1));
    F(i,:) = [fuelCost,activeLoss,emission,lvd];

    gPgLower = (cfg.constraint.pg_min-pg)./cfg.constraint.pg_scale;
    gPgUpper = (pg-cfg.constraint.pg_max)./cfg.constraint.pg_scale;
    gQgLower = (cfg.constraint.qg_min-qg)./cfg.constraint.qg_scale;
    gQgUpper = (qg-cfg.constraint.qg_max)./cfg.constraint.qg_scale;
    gVLower = (cfg.constraint.load_voltage_min-pqVoltage)./cfg.constraint.voltage_scale;
    gVUpper = (pqVoltage-cfg.constraint.load_voltage_max)./cfg.constraint.voltage_scale;
    G(i,:) = [gPgLower;gPgUpper;gQgLower;gQgUpper;gVLower;gVUpper]' - tolerance;

    detail.slack_pg(i)=result.gen(cfg.mapping.slack_gen_row,PG_COL);
    detail.total_pg(i)=sum(result.gen(result.gen(:,GEN_STATUS_COL)>0,PG_COL));
    detail.total_pd(i)=sum(result.bus(:,PD_COL));
    detail.minimum_pq_voltage(i)=min(pqVoltage);
    detail.maximum_pq_voltage(i)=max(pqVoltage);
    detail.cost_from_totcost(i)=fuelCost;
    detail.cost_from_coefficients(i)=costByCoeff;
    detail.cost_difference(i)=abs(fuelCost-costByCoeff);
end

positiveG=max(G,0);
detail.cv_sum=sum(positiveG,2);
detail.cv_max=max(positiveG,[],2);
detail.is_feasible=detail.success & detail.cv_max<=0;
detail.number_of_solutions=N;
detail.number_of_objectives=4;
detail.number_of_constraints=numberOfConstraints;
detail.objective_names=cfg.objective_names;
detail.decode_info=decodeInfo;
detail.constraint_sign_rule='G <= 0 feasible; G > 0 violated';
detail.constraint_normalization='groupwise physical range; tolerance subtracted';

cursor=0;
detail.constraint_index.pg_lower=cursor+(1:numberOfGenerators); cursor=cursor+numberOfGenerators;
detail.constraint_index.pg_upper=cursor+(1:numberOfGenerators); cursor=cursor+numberOfGenerators;
detail.constraint_index.qg_lower=cursor+(1:numberOfGenerators); cursor=cursor+numberOfGenerators;
detail.constraint_index.qg_upper=cursor+(1:numberOfGenerators); cursor=cursor+numberOfGenerators;
detail.constraint_index.voltage_lower=cursor+(1:numberOfPQ); cursor=cursor+numberOfPQ;
detail.constraint_index.voltage_upper=cursor+(1:numberOfPQ); cursor=cursor+numberOfPQ;
if cursor~=numberOfConstraints, error('IEEE legacy evaluation compatibility error.'); end
end


