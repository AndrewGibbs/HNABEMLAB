function [H,Hosc,HnonOsc] = besselhDecomp(NU,K,Z)
%similar to besselh.m, except decomposes Hankel function into oscillatory
%and non-oscillatory components, useful for oscillatory quadrature.

    if K ~=1
        error('Only first kind Hankel so far please');
    end
    
    threshSmallMed = -20;
    threshMedLarge = inf;
    smallInds = find(imag(Z)>threshSmallMed);
    medInds = find((imag(Z)<=threshSmallMed)) ;%.* (imag(Z)<=threshMedLarge));
    %largeInds = find(imag(Z)>=threshMedLarge);
    
    if ~isempty(smallInds)
        H(smallInds) = besselh(NU,K,Z(smallInds));
        Hosc(smallInds) = exp(1i*Z(smallInds));
        HnonOsc(smallInds) = H(smallInds)./Hosc(smallInds);
    end
    if ~isempty(medInds)
        [H(medInds), Hosc(medInds), HnonOsc(medInds)] = besselhExpn(NU,K,Z(medInds));
    end
%     if ~isempty(largeInds)
%         [H(largeInds), Hosc(largeInds), HnonOsc(largeInds)] = besselhAsmyp(NU, K, Z(largeInds));
%     end
    
    if ~isequal(size(Z),size(H))
        Zsize = size(Z);
        H = reshape(H,Zsize);
        Hosc = reshape(Hosc,Zsize);
        HnonOsc = reshape(HnonOsc,Zsize);
    end
end