function X=ieee118_fsa_decode(X,cfg)
    X=min(max(X,cfg.variable.lower),cfg.variable.upper);
    d=cfg.index.TAP;
    X(:,d)=cfg.variable.lower(d)+round((X(:,d)-cfg.variable.lower(d))/cfg.tapStep)*cfg.tapStep;
    X(:,d)=round(X(:,d),12);
    X(:,cfg.index.QC)=round(X(:,cfg.index.QC));
    X=min(max(X,cfg.variable.lower),cfg.variable.upper);
end
