function tf = opfpim3e_lexicographic_less(a,b)
%OPFPIU_LEXICOGRAPHIC_LESS Floating-point-aware lexicographic comparison.

    scale = max([ones(size(a));abs(a);abs(b)],[],1);
    tolerance = 16*eps(scale);
    tf = false;

    for k = 1:numel(a)
        if a(k) < b(k)-tolerance(k)
            tf = true;
            return;
        elseif a(k) > b(k)+tolerance(k)
            return;
        end
    end
end





