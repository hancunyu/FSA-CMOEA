function Dec = opfpim3e_project_dec(Problem,Dec,profile)
%OPFPIU_PROJECT_DEC Box projection and benchmark-defined discretization.

    if isempty(Dec), return; end
    Lower = repmat(Problem.lower,size(Dec,1),1);
    Upper = repmat(Problem.upper,size(Dec,1),1);
    Dec = min(max(Dec,Lower),Upper);

    if profile.mode == "ieee30_control"
        if ~isempty(profile.tapDims)
            d = profile.tapDims; step = profile.tapStep;
            Dec(:,d) = round(Dec(:,d)./step).*step;
        end
        if ~isempty(profile.integerDims)
            Dec(:,profile.integerDims) = round(Dec(:,profile.integerDims));
        end
        Dec = min(max(Dec,Lower),Upper);
    elseif profile.mode == "ieee57_control"
        if ~isempty(profile.tapDims)
            d = profile.tapDims; step = profile.tapStep;
            Dec(:,d) = Lower(:,d)+round((Dec(:,d)-Lower(:,d))./step).*step;
            Dec(:,d) = round(Dec(:,d),12);
        end
        if ~isempty(profile.integerDims)
            Dec(:,profile.integerDims) = round(Dec(:,profile.integerDims));
        end
        Dec = min(max(Dec,Lower),Upper);
    end
end





