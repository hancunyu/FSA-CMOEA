function cfg = ieee118_fsa_config()
% Deterministic extension of the manuscript IEEE57 control-space model.
% The other paper contributes audited network mappings, not RES objectives.
    persistent cached
    if ~isempty(cached), cfg=cached; return; end
    mpc=loadcase('case118');
    assert(size(mpc.bus,1)==118 && size(mpc.gen,1)==54 && size(mpc.branch,1)==186);
    gen=find(mpc.gen(:,8)>0);
    slackBus=mpc.bus(mpc.bus(:,2)==3,1);
    assert(isscalar(slackBus) && slackBus==69);
    pgRows=gen(mpc.gen(gen,1)~=slackBus);
    [ok,vgBusRows]=ismember(mpc.gen(gen,1),mpc.bus(:,1)); assert(all(ok));
    % Nine native non-unity ratios are adjustable; two unity-ratio entries
    % remain fixed, as in the audited mapping of the other IEEE118 model.
    tapRows=find(mpc.branch(:,11)>0 & ...
        ((mpc.branch(:,9)~=0 & abs(mpc.branch(:,9)-1)>eps) | mpc.branch(:,10)~=0));
    expectedTaps=[8;32;36;51;93;95;102;107;127];
    assert(isequal(tapRows,expectedTaps));
    qcBuses=[34;44;45;46;48;74;79;82;83;105;107;110];
    [ok,qcRows]=ismember(qcBuses,mpc.bus(:,1)); assert(all(ok));
    pqRows=find(mpc.bus(:,2)==1);
    assert(numel(pgRows)==53 && numel(pqRows)==64);
    cfg.baseMpc=mpc; cfg.D=128; cfg.tapStep=0.02;
    cfg.index.PG=1:53; cfg.index.VG=54:107;
    cfg.index.TAP=108:116; cfg.index.QC=117:128;
    cfg.mapping.pg_gen_rows=pgRows;
    cfg.mapping.vg_gen_rows=gen; cfg.mapping.vg_bus_rows=vgBusRows;
    cfg.mapping.tap_rows=tapRows; cfg.mapping.qc_bus_rows=qcRows;
    cfg.mapping.qc_buses=qcBuses; cfg.mapping.pq_bus_rows=pqRows;
    cfg.mapping.slack_gen_row=gen(mpc.gen(gen,1)==slackBus);
    cfg.variable.lower=[mpc.gen(pgRows,10);mpc.bus(vgBusRows,13);0.90*ones(9,1);zeros(12,1)]';
    cfg.variable.upper=[mpc.gen(pgRows,9);mpc.bus(vgBusRows,12);1.10*ones(9,1);20*ones(12,1)]';
    assert(all(cfg.variable.upper>cfg.variable.lower));
    cfg.constraint.pg_min=mpc.gen(gen,10);cfg.constraint.pg_max=mpc.gen(gen,9);
    cfg.constraint.qg_min=mpc.gen(gen,5); cfg.constraint.qg_max=mpc.gen(gen,4);
    cfg.constraint.v_min=mpc.bus(pqRows,13);cfg.constraint.v_max=mpc.bus(pqRows,12);
    cfg.constraint.pg_scale=max(cfg.constraint.pg_max-cfg.constraint.pg_min,1);
    cfg.constraint.qg_scale=max(cfg.constraint.qg_max-cfg.constraint.qg_min,1);
    cfg.constraint.v_scale=max(cfg.constraint.v_max-cfg.constraint.v_min,eps);
    cfg.constraint.tolerance=1e-6;cfg.constraint.penalty=1e6;
    cfg.constraint.count=344;
    cfg.constraint.groups={1:108,109:216,217:344};
    cfg.mpopt=mpoption('verbose',0,'out.all',0,'pf.alg','NR', ...
        'pf.tol',1e-8,'pf.nr.max_it',50,'pf.enforce_q_lims',0,'exp.use_legacy_core',1);
    cfg.activeBranches=find(mpc.branch(:,11)>0);
    assert(all(mpc.gencost(:,1)==2) && all(mpc.gencost(:,4)==3));
    cfg.cost=mpc.gencost(:,5:7);
    cached=cfg;
end
