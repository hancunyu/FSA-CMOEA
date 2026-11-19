function Offspring = opfpim3e_de_pbest1_dec(Problem,Parent1,Parent2,Parent3,Parent4,profile)
    P1=Parent1.decs; P2=Parent2.decs; P3=Parent3.decs; P4=Parent4.decs;
    [N,D]=size(P1);
    Fm=[0.6,0.8,1.0]; CRm=[0.1,0.2,1.0];
    F=Fm(randi(numel(Fm),N,1)); F=F'; F=F(:,ones(1,D));
    CR=CRm(randi(numel(CRm),N,1)); CR=CR'; CR=CR(:,ones(1,D));
    dims=opfpim3e_active_dims(Problem,profile,false);
    Site=false(N,D); Site(:,dims)=rand(N,numel(dims))<CR(:,dims);
    Offspring=P1;
    base  = Offspring(Site);
    scale = F(Site);
    delta = P2(Site)-Offspring(Site)+P3(Site)-P4(Site);
    Offspring(Site)=base(:)+scale(:).*delta(:);
    Offspring=opfpim3e_polynomial_mutation(Problem,Offspring,dims);
end





