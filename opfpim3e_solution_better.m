function tf = opfpim3e_solution_better(a,b,profile,objLo,objScale)
    fa = all(a.cons<=0,2);
    fb = all(b.cons<=0,2);

    if fa && ~fb
        tf = true;
        return;
    elseif ~fa && fb
        tf = false;
        return;
    end

    ka = opfpim3e_solution_key(a,profile,objLo,objScale);
    kb = opfpim3e_solution_key(b,profile,objLo,objScale);
    tf = opfpim3e_lexicographic_less(ka,kb);
end





