function profile = opfpim3e_validate_problem(Problem)
%OPFPIU_VALIDATE_PROBLEM Validate and classify the OPF problem interface.

    if isempty(Problem.D) || Problem.D<1 || Problem.D~=round(Problem.D)
        error('OPFPI_M3E:InvalidDimension','Problem.D must be a positive integer.');
    end
    if numel(Problem.lower)~=Problem.D || numel(Problem.upper)~=Problem.D
        error('OPFPI_M3E:InvalidBounds','The number of lower/upper bounds must equal Problem.D.');
    end
    if any(~isfinite(Problem.lower)) || any(~isfinite(Problem.upper)) || any(Problem.upper<Problem.lower)
        error('OPFPI_M3E:InvalidBounds','Decision bounds are invalid.');
    end

    profile = opfpim3e_problem_profile(Problem);
    if profile.mode=="legacy_rwmop" && Problem.D~=34
        error('OPFPI_M3E:LegacyDimension','The frozen RWMOP restoration requires D=34.');
    end
    if profile.mode=="ieee30_control" && Problem.D~=24
        error('OPFPI_M3E:IEEE30Dimension','IEEE30_MOOPF control encoding requires D=24.');
    end
    if profile.mode=="ieee57_control" && Problem.D~=33
        error('OPFPI_M3E:IEEE57Dimension','IEEE57_MOOPF control encoding requires D=33.');
    end
end





