function [return_pop,return_Fitness] = opfpim3e_aux_selection_original(Population,N,VAR)
% Semantics aligned with the supplied IMTCMO auxiliary selection.
    cv=sum(max(0,Population.cons),2);
    relaxed=cv<=VAR;
    fPopulation=Population(relaxed); ifPopulation=Population(~relaxed);
    fFitness=[]; ifFitness=[];
    if isempty(fPopulation)
        [ifPopulation,ifFitness]=opfpim3e_main_selection(ifPopulation,N,true);
    elseif length(fPopulation)<=N
        aug=[fPopulation.objs,cv(relaxed)];
        fit=opfpim3e_calfitness(aug);
        [~,ord]=sort(fit); fPopulation=fPopulation(ord); fFitness=fit(ord);
        need=N-length(fPopulation);
        if need>0 && ~isempty(ifPopulation)
            [ifPopulation,ifFitness]=opfpim3e_main_selection(ifPopulation,min(need,length(ifPopulation)),true);
            if ~isempty(fFitness), ifFitness=ifFitness+max(fFitness); end
        else
            ifPopulation=[]; ifFitness=[];
        end
    else
        aug=[fPopulation.objs,cv(relaxed)];
        fit=opfpim3e_calfitness(aug);
        Next=fit<1;
        if sum(Next)<N
            [~,ord]=sort(fit); Next(ord(1:N))=true;
        elseif sum(Next)>N
            Del=opfpim3e_truncation(fPopulation(Next).objs,sum(Next)-N);
            tmp=find(Next); Next(tmp(Del))=false;
        end
        fPopulation=fPopulation(Next); fFitness=fit(Next);
        [fFitness,ord]=sort(fFitness); fPopulation=fPopulation(ord);
        ifPopulation=[]; ifFitness=[];
    end
    return_pop=[fPopulation,ifPopulation];
    return_Fitness=[fFitness,ifFitness];
end





