function [ca, cb, cxi] = hoi_asy_I_sp1(s, a, b, xi, fa, fb, fxi, ga, gb, gxi)
%function [ca, cb, cxi] = hoi_asy_I_sp1(s, a, b, xi, fa, fb, fxi, ga, gb, gxi)
%
%   Return the asymptotic expansion of I in the presence of a stationary
%   point of order 1.
%   The method of asymptotic order s requires the derivatives up to order
%   2s-2 for f, and up to order 2s-1 for g.

% $Id: hoi_asy_I_sp1.m 2 2006-11-03 09:18:49Z daan $

ca = zeros(s,1);
cb = zeros(s,1);
cxi = zeros(s,1);

rho_xi = hoi_asy_rho_xi(fxi, gxi, s-1);
rho_a = hoi_asy_rho(fa, ga, rho_xi, s-1);
rho_b = hoi_asy_rho(fb, gb, rho_xi, s-1);

for m=1:s
    cxi(m) = (1i)^(m-1)*rho_xi(m);
    ca(m) = -(1i)^m*(rho_a(m)-rho_xi(m))/ga(2);
    cb(m) = -(1i)^m*(rho_b(m)-rho_xi(m))/gb(2);
end
