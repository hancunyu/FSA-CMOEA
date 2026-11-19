function [H,J,meta] = opfpim3e_legacy_signed_residual_jacobian(Dec,needJacobian)
%OPFPI_C2_SIGNED_RESIDUAL_JACOBIAN Signed OPF mismatch and analytic state Jacobian.
%
% Decision layout:
%   1:13  = real(V2:V14)
%   14:26 = imag(V2:V14)
%   27:30 = Pg at buses [2,3,6,8]
%   31:34 = Qg at buses [2,3,6,8]
%
% H contains [delP(2:14),delQ(2:14)] before abs(.)-1e-4.
% J is dH/d[real(V2:V14),imag(V2:V14)] and has size 26x26xN.
%
% This function does not call Problem.Evaluation and therefore consumes no
% FE. All residual and Jacobian calls are separately counted by OPFPI_C2.

    if nargin < 2
        needJacobian = nargout >= 2;
    end
    if isempty(Dec)
        H = zeros(0,26);
        J = zeros(26,26,0);
        meta = struct('voltageMagnitude',zeros(0,13));
        return;
    end
    if size(Dec,2) ~= 34
        error('OPFPI_C2:InvalidDecisionDimension',...
            'Expected a decision matrix with 34 columns.');
    end

    [Y,Pload,Qload,G,B] = opfpim3e_legacy_pf14_data();
    N = size(Dec,1);
    H = zeros(N,26);
    if needJacobian
        J = zeros(26,26,N);
    else
        J = [];
    end
    voltageMagnitude = zeros(N,13);
    genBus = [2,3,6,8];

    for n = 1:N
        e = [1;Dec(n,1:13)'];
        f = [0;Dec(n,14:26)'];
        V = e + 1j*f;

        Pg = zeros(14,1);
        Qg = zeros(14,1);
        Pg(genBus) = Dec(n,27:30)';
        Qg(genBus) = Dec(n,31:34)';

        I = Y*V;
        a = real(I);
        b = imag(I);
        S = V.*conj(I);
        Psp = real(S);
        Qsp = imag(S);

        delP = Psp - Pg + Pload;
        delQ = Qsp - Qg + Qload;
        H(n,:) = [delP(2:14)',delQ(2:14)'];
        voltageMagnitude(n,:) = abs(V(2:14))';

        if needJacobian
            Jn = zeros(26,26);
            for ii = 1:13
                i = ii + 1;
                for kk = 1:13
                    k = kk + 1;
                    dik = double(i == k);

                    % P_i = e_i*a_i + f_i*b_i
                    Jn(ii,kk) = dik*a(i) + e(i)*G(i,k) + f(i)*B(i,k);
                    Jn(ii,13+kk) = dik*b(i) - e(i)*B(i,k) + f(i)*G(i,k);

                    % Q_i = f_i*a_i - e_i*b_i
                    Jn(13+ii,kk) = f(i)*G(i,k) - dik*b(i) - e(i)*B(i,k);
                    Jn(13+ii,13+kk) = dik*a(i) - f(i)*B(i,k) - e(i)*G(i,k);
                end
            end
            J(:,:,n) = Jn;
        end
    end

    meta = struct('voltageMagnitude',voltageMagnitude);
end





