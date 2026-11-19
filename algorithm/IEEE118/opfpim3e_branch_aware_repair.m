function [RepairedMain,RepairedAux,diag] = opfpim3e_branch_aware_repair(...
    Problem,Samples,Candidates,profile)
%opfpim3e_BRANCH_AWARE_REPAIR Dispatch restoration by representation.
%
% For IEEE30/57 control-space MOOPF, ODP returns two repaired archives:
%   RepairedMain: repairs that preserve objective direction and may enter
%                 the strict feasibility-first main population.
%   RepairedAux : all accepted repairs, including objective-direction-
%                 damaged repairs, for the relaxed auxiliary population.
%
% This soft demotion changes no objective values and adds no extra FE.

    if profile.mode == "legacy_rwmop"
        [R,diag] = opfpim3e_legacy_branch_repair(Problem,Candidates);
        RepairedMain = R;
        RepairedAux  = R;
    else
        [RepairedMain,RepairedAux,diag] = opfpim3e_control_branch_repair(...
            Problem,Samples,Candidates,profile);
    end
end





