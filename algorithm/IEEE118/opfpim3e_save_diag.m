function opfpim3e_save_diag(Algorithm,Problem,diag)
% Diagnostic storage only; no search, evaluation or RNG effects.
    Algorithm.diagnostics = diag;
end
