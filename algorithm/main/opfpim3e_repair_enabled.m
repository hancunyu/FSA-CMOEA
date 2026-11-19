function enabled = opfpim3e_repair_enabled(Population1,profile)
%OPFPIU_REPAIR_ENABLED Representation-aware physical-restoration gate.
% RWMOP restoration remains active. For IEEE control-only encodings,
% restoration is disabled once the published main population is all feasible.

    controlModes = ["ieee30_control","ieee57_control"];
    if ~ismember(profile.mode,controlModes)
        enabled = true;
        return;
    end
    if isempty(Population1) || isempty(Population1.cons)
        enabled = false;
        return;
    end
    enabled = ~all(all(Population1.cons<=0,2));
end





