function diag = opfpim3e_diag_snapshot(diag,Problem,P1,P2,Published)
    C = [P1.cons;P2.cons];
    V = max(0,C);
    k = numel(diag.trace.FE)+1;

    diag.trace.FE(k,1) = Problem.FE;
    diag.trace.feasible1(k,1) = sum(all(P1.cons<=0,2));
    diag.trace.feasible2(k,1) = sum(all(P2.cons<=0,2));
    diag.trace.publishedFeasible(k,1) = ...
        sum(all(Published.cons<=0,2));

    if isempty(V)
        diag.trace.bestCVsum(k,1) = NaN;
        diag.trace.bestCVmax(k,1) = NaN;
    else
        diag.trace.bestCVsum(k,1) = min(sum(V,2));
        diag.trace.bestCVmax(k,1) = min(max(V,[],2));
    end

    diag.trace.population1SHA{k,1} = opfpim3e_sha256_matrix(P1.decs);
    diag.trace.population2SHA{k,1} = opfpim3e_sha256_matrix(P2.decs);
    diag.finalFE = Problem.FE;
end





