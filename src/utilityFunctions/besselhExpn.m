function [full,expBit,sum] = besselhExpn(v, kind, z)
%computes asymptotic expansion of Hankel function of first kind order v
% 
%     if nargin == 2
%         kind = 1;
%     end
    
    %doesn't work so well for negative phase, so use DLMF analytic
    %continuation identity:
    if min(min(angle(z)))<0
        %split into negtive and positive bits
        [full,expBit,sum] = besselhAsmypConj(@(K,w) besselhExpn(v, K, w), kind, z);
        return;
    end
    
    pm = (-1)^(kind+1);
    
    thresh = 1e-16;
    summand = inf;
    k=0;
    
    expBit = exp(pm*1i*z);
    preFactor = sqrt(2./(pi*z)) * exp(pm*1i*(- .5*pi*v - .25*pi));
    sum = 0;
    
    keepAddingTerms = true;
    while keepAddingTerms %min(min(abs(summand)))>thresh
        summandPrev = summand;
        summand = (pm*1i)^k * a_k(k) * z.^-k;
        sum = sum + summand;
        if k == 0
            thresh = abs(max(sum))*1e-16;
        end
        if abs(max(summand)) > abs(max(summandPrev)) || abs(max(summand))<thresh
            keepAddingTerms = false;
        end
        k = k+1;
    end
 
    full = expBit .* preFactor .* sum;
    
    function A = a_k(k)
        if k == 0
            A = 1;
        else
            A = 1/(factorial(k)*8^k);
            for j = 1:2:(2*k-1)
                A = A *(4*v^2 - j^2);
            end
        end
    end

end
