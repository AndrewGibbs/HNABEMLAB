function vals = besselh_0_1_nonosc_large_imag(z)
% Matlab's besselh(0,1,z) can become inaccurate when Im(z) is large.
% This is particularly problematic when dividing through by exp(i*z).
% This code deals with the large Im(z) entries using homemade code.

    Im_thresh = 30; %based on quick experiments
    z_matlab_inds = imag(z)<Im_thresh;
    z_asymp_inds = imag(z)>=Im_thresh;
    vals = zeros(size(z));
    vals(z_matlab_inds) = besselh(0,1,z(z_matlab_inds))./exp(1i*z(z_matlab_inds));
    [~, ~, HnonOsc] = besselhAsmyp(0, 1, z(z_asymp_inds));
    vals(z_asymp_inds) = HnonOsc;
end