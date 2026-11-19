function DecOut = opfpim3e_legacy_flat_start_dec(DecIn,Problem)
%OPFPI_C3_FLAT_START_DEC Deterministic high-voltage flat-state restart.
%
% State:
%   Re(V2:V14) = 1
%   Im(V2:V14) = 0
%
% Controls 27:34 are copied exactly from the input candidate.
% The decision is clipped only to the benchmark box bounds.

    if isempty(DecIn)
        DecOut = DecIn;
        return;
    end
    if size(DecIn,2) ~= 34
        error('OPFPI_C3:InvalidDecisionDimension',...
            'Expected 34 decision variables.');
    end

    DecOut = DecIn;
    DecOut(:,1:13) = 1;
    DecOut(:,14:26) = 0;
    DecOut = max(min(DecOut,Problem.upper),Problem.lower);
end





