function [Population,Fitness] = opfpim3e_main_selection(Population,N,isOrigin)
    if isOrigin
        Fitness=opfpim3e_calfitness(Population.objs,Population.cons);
    else
        Fitness=opfpim3e_calfitness(Population.objs);
    end
    Next=Fitness<1;
    if sum(Next)<N
        [~,Rank]=sort(Fitness); Next(Rank(1:min(N,numel(Rank))))=true;
    elseif sum(Next)>N
        Del=opfpim3e_truncation(Population(Next).objs,sum(Next)-N);
        Temp=find(Next); Next(Temp(Del))=false;
    end
    Population=Population(Next); Fitness=Fitness(Next);
    [Fitness,rank]=sort(Fitness); Population=Population(rank);
end





