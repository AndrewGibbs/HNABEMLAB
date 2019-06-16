function [ca,cb] = hoi_exp_I(s, a, b, fa, fb, ga, gb)
%function [ca,cb] = hoi_exp_I(s, a, b, fa, fb, ga, gb)
%
%   Return the full asymptotic expansion of I in the absence of stationary
%   points.
%   The method of asymptotic order s requires the derivatives up to order
%   s-1 for f, and up to order s for g.

% $Id: hoi_exp_I.m 2 2006-11-03 09:18:49Z daan $

ca = hoi_exp_F(s, a, fa, ga);
cb = hoi_exp_F(s, b, fb, gb);
