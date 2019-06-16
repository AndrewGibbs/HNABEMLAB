function [ca,cb,cxi] = hoi_exp_I_sp1(s, n, a, b, xi, fa, fb, fxi, ga, gb, gxi)
%function [ca,cb,cxi] = hoi_exp_I_sp1(s, n, a, b, xi, fa, fb, fxi, ga, gb, gxi)
%
%   Return the full asymptotic expansion of I in the presence of a
%   stationary point of order 1 at xi.
%   The method of asymptotic order s requires the derivatives up to order
%   s-1 for f, and up to order s for g. The expansion around xi requires n
%   derivatives of f, and n+1 derivatives of g.
%   A match of orders requires that s+1 = (n+1)/2.

% $Id: hoi_exp_I_sp1.m 2 2006-11-03 09:18:49Z daan $

ca = hoi_exp_F(s, a, fa, ga);
cb = hoi_exp_F(s, b, fb, gb);
[c1,c2] = hoi_exp_F_sp1(n, xi, fxi, gxi);
cxi = c1-c2;
