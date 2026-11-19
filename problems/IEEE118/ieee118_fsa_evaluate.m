function [F,G,detail,X]=ieee118_fsa_evaluate(X,cfg)
    if nargin<2,cfg=ieee118_fsa_config();end
    X=ieee118_fsa_decode(X,cfg);
    n=size(X,1);F=ones(n,2)*cfg.constraint.penalty;
    G=ones(n,cfg.constraint.count)*cfg.constraint.penalty;
    detail.success=false(n,1);detail.errors=strings(n,1);
    for i=1:n
        mpc=cfg.baseMpc;x=X(i,:);
        mpc.gen(cfg.mapping.pg_gen_rows,2)=x(cfg.index.PG)';
        mpc.gen(cfg.mapping.vg_gen_rows,6)=x(cfg.index.VG)';
        mpc.bus(cfg.mapping.vg_bus_rows,8)=x(cfg.index.VG)';
        mpc.branch(cfg.mapping.tap_rows,9)=x(cfg.index.TAP)';
        mpc.bus(cfg.mapping.qc_bus_rows,6)=x(cfg.index.QC)';
        try
            [r,success]=runpf(mpc,cfg.mpopt);
        catch err
            detail.errors(i)=string(err.message);continue;
        end
        if ~success || any(~isfinite(r.bus(:))) || any(~isfinite(r.gen(:)))
            continue;
        end
        pg=r.gen(cfg.mapping.vg_gen_rows,2);
        qg=r.gen(cfg.mapping.vg_gen_rows,3);
        v=r.bus(cfg.mapping.pq_bus_rows,8);
        c=cfg.cost;
        cost=sum(c(:,1).*r.gen(:,2).^2+c(:,2).*r.gen(:,2)+c(:,3));
        loss=sum(r.branch(cfg.activeBranches,14)+r.branch(cfg.activeBranches,16));
        if loss<0 && loss>-1e-8,loss=0;end
        gi=[(cfg.constraint.pg_min-pg)./cfg.constraint.pg_scale;
            (pg-cfg.constraint.pg_max)./cfg.constraint.pg_scale;
            (cfg.constraint.qg_min-qg)./cfg.constraint.qg_scale;
            (qg-cfg.constraint.qg_max)./cfg.constraint.qg_scale;
            (cfg.constraint.v_min-v)./cfg.constraint.v_scale;
            (v-cfg.constraint.v_max)./cfg.constraint.v_scale]'-cfg.constraint.tolerance;
        if any(~isfinite([cost,loss,gi])),continue;end
        F(i,:)=[cost,loss];G(i,:)=gi;detail.success(i)=true;
    end
    detail.feasible=detail.success & all(G<=0,2);
end
