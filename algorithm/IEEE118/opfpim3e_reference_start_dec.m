function DecOut = opfpim3e_reference_start_dec(DecIn,Problem,profile)
%OPFPIU_REFERENCE_START_DEC Deterministic physical reference path.
% PG is preserved for IEEE control-only problems; voltage/tap controls are
% reset toward 1 p.u. and shunt compensation toward 0.

    DecOut = DecIn;
    if isempty(DecIn), return; end

    if profile.mode == "ieee30_control"
        DecOut(:,6:11) = 1;
        DecOut(:,12:15) = 1;
        DecOut(:,16:24) = 0;
    elseif profile.mode == "ieee57_control"
        DecOut(:,7:13) = 1;
        DecOut(:,14:30) = 1;
        DecOut(:,31:33) = 0;
    elseif profile.mode == "ieee118_control"
        cfg = Problem.Config;
        DecOut(:,cfg.index.VG) = 1;
        DecOut(:,cfg.index.TAP) = 1;
        DecOut(:,cfg.index.QC) = 0;
    else
        midpoint = (Problem.lower+Problem.upper)/2;
        DecOut(:,profile.referenceDims) = repmat(midpoint(profile.referenceDims),size(DecIn,1),1);
    end

    DecOut = opfpim3e_project_dec(Problem,DecOut,profile);
end





