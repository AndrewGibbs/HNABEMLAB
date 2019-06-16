function [ca, cb] = hoi_asy_I(s, a, b, fa, fb, ga, gb)
%function [ca, cb] = hoi_asy_I(s, a, b, fa, fb, ga, gb)
%
%   Return the asymptotic expansion of I in the absence of stationary
%   points.
%   The method of asymptotic order s requires the derivatives up to order
%   s-1 for f, and up to order s for g.

% $Id: hoi_asy_I.m 2 2006-11-03 09:18:49Z daan $

ca = zeros(s,1);
cb = zeros(s,1);

sigma_a = hoi_asy_sigma(fa, ga, s-1);
sigma_b = hoi_asy_sigma(fb, gb, s-1);

for m=1:s
    ca(m) = -(1i)^m*sigma_a(m)/ga(2);
    cb(m) = -(1i)^m*sigma_b(m)/gb(2);
end
