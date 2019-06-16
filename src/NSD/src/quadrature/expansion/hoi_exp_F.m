function c = hoi_exp_F(s, a, fa, ga)
%function c = hoi_exp_F(s, a, fa, ga)
%
%   Return the full asymptotic expansion of F at a regular point a.
%   The expansion of asymptotic order s requires the derivatives up to order
%   s-1 for f, and up to order s for g.

% $Id: hoi_exp_F.m 2 2006-11-03 09:18:49Z daan $

c = zeros(s,1);

sigma_a = hoi_asy_sigma(fa, ga, s-1);

for m=1:s
    c(m) = -(1i)^m*sigma_a(m)/ga(2);
end
