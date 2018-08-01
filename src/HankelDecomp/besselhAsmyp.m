function [H, Hosc, HnonOsc] = besselhAsmyp(v, kind, z)
%a simpler asymptotic approximation

    if min(min(angle(z)))<0 && kind == 1
        %split into negtive and positive bits
        [H, Hosc, HnonOsc] = besselhAsmypConj(@(K,w) besselhAsmyp(v, K, w), kind, z);
        return;
    end
    pm = (-1)^(kind+1);
    Hosc = exp(pm*1i*z);
    HnonOsc = sqrt(2./(pi*z)) * exp(pm*1i*(- v*pi/2 - pi/4));
    H = HnonOsc .* Hosc;
end

