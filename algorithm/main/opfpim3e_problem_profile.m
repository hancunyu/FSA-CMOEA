function profile = opfpim3e_problem_profile(Problem)
%OPFPIU_PROBLEM_PROFILE Problem-specific data only; no search parameters.

    name = string(class(Problem));
    profile = struct();
    profile.problemClass = name;
    profile.mode = "generic_control";
    profile.controlDims = 1:Problem.D;
    profile.preserveDims = [];
    profile.referenceDims = 1:Problem.D;
    profile.constraintGroups = {};
    profile.constraintGroupNames = strings(0,1);
    profile.tapDims = [];
    profile.tapStep = [];
    profile.integerDims = [];

    if Problem.D == 34 && ~isempty(regexp(char(name),'^RWMOP4[0-6]$','once'))
        profile.mode = "legacy_rwmop";
        profile.controlDims = 27:34;
        profile.preserveDims = 27:34;
        profile.referenceDims = 1:26;
        return;
    end

    if Problem.D == 24 && contains(upper(name),"IEEE30")
        profile.mode = "ieee30_control";
        profile.controlDims = 1:24;
        profile.preserveDims = 1:5;
        profile.referenceDims = 6:24;
        profile.tapDims = 12:15;
        profile.tapStep = 0.02;
        profile.integerDims = 16:24;
        profile.constraintGroups = {1:12,13:24,25:72,73:113};
        profile.constraintGroupNames = ["PG";"QG";"Voltage";"Branch"];
        return;
    end

    if Problem.D == 33 && contains(upper(name),"IEEE57")
        profile.mode = "ieee57_control";
        profile.controlDims = 1:33;
        profile.preserveDims = 1:6;
        profile.referenceDims = 7:33;
        profile.tapDims = 14:30;
        profile.tapStep = 0.02;
        profile.integerDims = 31:33;
        profile.constraintGroups = {1:14,15:28,29:128};
        profile.constraintGroupNames = ["PG";"QG";"Voltage"];
        return;
    end
end





