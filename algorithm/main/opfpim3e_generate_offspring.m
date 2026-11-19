function Offspring = opfpim3e_generate_offspring(Problem,Population,...
    Fitness,Zmin,nOff,profile)
% Generate up to nOff candidates with the frozen C3 variation mixture.

    Offspring = [];
    if nOff <= 0
        return;
    end

    popCount = length(Population);
    order = randperm(popCount);
    MatingPool = Population(order);

    [Mate1,Mate2,Mate3] = opfpim3e_neighbor_pairing(MatingPool,Zmin);

    nRand = floor(nOff/2);
    nPbest = nOff-nRand;
    Dec = [];

    if nRand > 0
        Dec = [Dec;opfpim3e_de_rand1_dec(Problem,...
            Mate1(1:nRand),Mate2(1:nRand),Mate3(1:nRand),profile)]; %#ok<AGROW>
    end
    if nPbest > 0
        Dec = [Dec;opfpim3e_de_pbest1_main_dec(Population,nPbest,...
            Problem,Fitness,0.1,profile)]; %#ok<AGROW>
    end

    Dec = opfpim3e_project_dec(Problem,Dec,profile);
    Offspring = Problem.Evaluation(Dec);
end





