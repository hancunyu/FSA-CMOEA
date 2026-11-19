function diag = opfpim3e_merge_diag(diag,increment)
%OPFPIU_MERGE_DIAG Merge numeric diagnostic counters safely.

    if isempty(increment)
        return;
    end
    names = fieldnames(increment);

    for i = 1:numel(names)
        name = names{i};
        value = increment.(name);

        if strcmp(name,'firstRepairCallFE')
            if (~isfield(diag,name) || isnan(diag.(name))) && ...
                    isfinite(value)
                diag.(name) = value;
            end
            continue;
        end

        if strcmp(name,'conditionNumberMax')
            if ~isfield(diag,name)
                diag.(name) = value;
            else
                diag.(name) = max(diag.(name),value);
            end
            continue;
        end

        if isnumeric(value) && isscalar(value)
            if ~isfield(diag,name)
                diag.(name) = value;
            else
                diag.(name) = diag.(name)+value;
            end
        else
            diag.(name) = value;
        end
    end
end





