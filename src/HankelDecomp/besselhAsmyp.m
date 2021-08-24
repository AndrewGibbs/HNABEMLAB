function [H, Hosc, HnonOsc] = besselhAsmyp(v, kind, z)
if isempty(z)
    H=[]; Hosc=[]; HnonOsc=[];
    return;
end
thresh = eps;
kmax = 16;%ceil(log(2/pi/thresh)/log(max(max(abs(z))))-1/2);
    if min(min(angle(z)))<0 && kind == 1
        %split into negtive and positive bits
        [H, Hosc, HnonOsc] = besselhAsmypConj(@(K,w) besselhAsmyp(v, K, w), kind, z);
        return;
    end
    pm = (-1)^(kind+1);
    Hosc = exp(pm*1i*z);    
    % proceed using https://dlmf.nist.gov/10.17.
    % the below loop is designed to continue summing each for each entry z, 
    % in such a way where the sum stops accumulating when it begins to
    % grow, noting the divergence of the asymptotic expansion.
    
    Asum = ones(size(z));
    sum_old = Asum;
    keep_inds = true(size(z));
    sum_new = zeros(size(z));
    for k=1:kmax
        sum_new(keep_inds) = 1i^k*Hankel_asym_coeffs(k,v)./z(keep_inds).^k;
        %only continue summing at indices which are below desired
        %threshold, and for which the summand is still decreasing.
        keep_inds = (abs(sum_new)<abs(sum_old) | abs(sum_new)>thresh) & keep_inds;
        if sum(keep_inds)==0
            break;
        end
        Asum(keep_inds) = Asum(keep_inds) + sum_new(keep_inds);
        sum_old = sum_new;
    end
    HnonOsc = sqrt(2./(pi*z)) .* exp(pm*1i*(- v*pi/2 - pi/4)) .* Asum;
    H = HnonOsc .* Hosc;
end

