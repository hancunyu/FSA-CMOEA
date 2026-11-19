function dims = opfpim3e_active_dims(Problem,profile,controlOnly)
    if nargin < 3 || ~controlOnly
        dims = 1:Problem.D;
    else
        dims = profile.controlDims;
    end
end





