function [ca, cxi] = hoi_exp_I_sp1b(s, n, a, xi, fa, fxi, ga, gxi)
%function [ca, cxi] = hoi_exp_I_sp1b(s, n, a, xi, fa, fxi, ga, gxi)
%
%   Return the full asymptotic expansion of I in the presence of a
%   stationary point of order 1 at the endpoint b.
%   The method of asymptotic order s requires the derivatives up to order
%   s-1 for f, and up to order s for g. The expansion around xi requires n
%   derivatives of f, and n+1 derivatives of g.

% $Id: hoi_exp_I_sp1b.m 2 2006-11-03 09:18:49Z daan $

ca = hoi_exp_F(s, a, fa, ga);
[cxi,dummy] = hoi_exp_F_sp1(n, xi, fxi, gxi);
