function [r1,r2,r3] = opfpim3e_gn_indices(NP1,r0)
    NP0 = length(r0);
    r1 = floor(rand(1,NP0)*NP1)+1;
    for k = 1:999
        pos = r1==r0; if ~any(pos), break; end
        r1(pos)=floor(rand(1,sum(pos))*NP1)+1;
    end
    r2 = floor(rand(1,NP0)*NP1)+1;
    for k = 1:999
        pos = r2==r1 | r2==r0; if ~any(pos), break; end
        r2(pos)=floor(rand(1,sum(pos))*NP1)+1;
    end
    r3 = floor(rand(1,NP0)*NP1)+1;
    for k = 1:999
        pos = r3==r1 | r3==r0 | r3==r2; if ~any(pos), break; end
        r3(pos)=floor(rand(1,sum(pos))*NP1)+1;
    end
end





