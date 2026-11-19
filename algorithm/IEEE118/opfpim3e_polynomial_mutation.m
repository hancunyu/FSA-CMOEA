function Offspring = opfpim3e_polynomial_mutation(Problem,Offspring,dims)
    [N,D]=size(Offspring); disM=20;
    Lower=repmat(Problem.lower,N,1); Upper=repmat(Problem.upper,N,1);
    Offspring=min(max(Offspring,Lower),Upper);
    Site=false(N,D);
    Site(:,dims)=rand(N,numel(dims)) < 1/max(1,numel(dims));
    mu=rand(N,D);
    temp=Site & mu<=0.5;
    den=Upper-Lower; den(den==0)=1;
    Offspring(temp)=Offspring(temp)+den(temp).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(Offspring(temp)-Lower(temp))./den(temp)).^(disM+1)).^(1/(disM+1))-1);
    temp=Site & mu>0.5;
    Offspring(temp)=Offspring(temp)+den(temp).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(Upper(temp)-Offspring(temp))./den(temp)).^(disM+1)).^(1/(disM+1)));
    Offspring=min(max(Offspring,Lower),Upper);
end





