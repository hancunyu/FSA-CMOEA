function key = opfpim3e_legacy_branch_key(Dec)
%OPFPI_C3_BRANCH_KEY Branch-quality key for strictly feasible solutions.
%
% Lower is better:
%   1. voltage deviation from 1 p.u.
%   2. negative minimum voltage magnitude
%   3. maximum absolute voltage deviation
%
% This key is used only after strict feasibility has been achieved. It does
% not modify the RWMOP constraints or objective vector.

    if isempty(Dec)
        key = [Inf,Inf,Inf];
        return;
    end
    e = Dec(1:13);
    f = Dec(14:26);
    vm = sqrt(e.^2+f.^2);
    key = [sum((1-vm).^2),-min(vm),max(abs(1-vm))];
end





