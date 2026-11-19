function opfpim3e_save_diag(Algorithm,Problem,diag)
%opfpim3e_SAVE_DIAG Save diagnostics without affecting the search.

    folder = fullfile(pwd,'Data','OPFPI_Diagnostics');
    if ~exist(folder,'dir')
        mkdir(folder);
    end

    run = Algorithm.run;
    if isempty(run)
        run = 0;
    end

    stamp = char(javaMethod('randomUUID','java.util.UUID'));
    file = fullfile(folder,sprintf(...
        'OPFPI_M3E_%s_run%03d_%s_diag.mat',class(Problem),run,stamp));
    tmp = [tempname(folder),'.mat'];
    try
        save(tmp,'diag','-v7');
        movefile(tmp,file,'f');
    catch
        if exist(tmp,'file')
            delete(tmp);
        end
    end
end





