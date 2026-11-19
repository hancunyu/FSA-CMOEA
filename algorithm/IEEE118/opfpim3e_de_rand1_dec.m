function Offspring = opfpim3e_de_rand1_dec(Problem,Parent1,Parent2,Parent3,profile)
% Return decisions only; evaluation is explicit in OPFPI for FE accounting.
    Parent1 = Parent1.decs; Parent2 = Parent2.decs; Parent3 = Parent3.decs;
    [N,D] = size(Parent1);
    Fm = [0.6,0.8,1.0]; CRm = [0.1,0.2,1.0];
    F = Fm(randi(numel(Fm),N,1)); F = F'; F = F(:,ones(1,D));
    CR = CRm(randi(numel(CRm),N,1)); CR = CR'; CR = CR(:,ones(1,D));
    dims = opfpim3e_active_dims(Problem,profile,false);
    Site = false(N,D); Site(:,dims) = rand(N,numel(dims)) < CR(:,dims);
    Offspring = Parent1;
    base  = Offspring(Site);
    scale = F(Site);
    delta = Parent2(Site)-Parent3(Site);
    Offspring(Site) = base(:)+scale(:).*delta(:);
    Offspring = opfpim3e_polynomial_mutation(Problem,Offspring,dims);
end





