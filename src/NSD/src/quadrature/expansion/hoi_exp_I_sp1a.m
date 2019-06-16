function [cxi,cb] = hoi_exp_I_sp1a(s, n, xi, b, fxi, fb, gxi, gb)
%function [cxi,cb] = hoi_exp_I_sp1a(s, n, xi, b, fxi, fb, gxi, gb)
%
%   Return the full asymptotic expansion of I in the presence of a
%   stationary point of order 1 at the endpoint a.
%   The method of asymptotic order s requires the derivatives up to order
%   s-1 for f, and up to order s for g. The expansion around xi requires n
%   derivatives of f, and n+1 derivatives of g.

% $Id: hoi_exp_I_sp1a.m 2 2006-11-03 09:18:49Z daan $

[c1,cxi] = hoi_exp_F_sp1(n, xi, fxi, gxi);
cb = hoi_exp_F(s, b, fb, gb);
