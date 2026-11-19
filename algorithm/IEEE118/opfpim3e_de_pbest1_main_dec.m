function Offspring = opfpim3e_de_pbest1_main_dec(Population,popsize,Problem,Fitness,p,profile)
    permutation = randperm(Problem.N);
    r0 = permutation;
    [r1,r2,~] = opfpim3e_gn_indices(Problem.N,r0);
    array = permutation(1:popsize);
    [~,indBest] = sort(Fitness,'ascend');
    pNP = max(round(p*Problem.N),2);
    randindex = max(1,ceil(rand(1,popsize)*pNP));
    pbest = Population(indBest(randindex));
    Offspring = opfpim3e_de_pbest1_dec(Problem,Population(array),pbest,...
        Population(r1(1:popsize)),Population(r2(1:popsize)),profile);
end





