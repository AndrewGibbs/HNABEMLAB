function [ca, cxi1, cxi2] = hoi_asy_I_sp1b(s, a, xi, fa, fxi, ga, gxi)
%function [ca, cxi1, cxi2] = hoi_asy_I_sp1b(s, a, xi, fa, fxi, ga, gxi)
%
%   Return the asymptotic expansion of I in the presence of a stationary
%   point of order 1.
%   The method of asymptotic order s requires the derivatives up to order
%   2s-1 for f, and up to order 2s for g.

% $Id: hoi_asy_I_sp1b.m 2 2006-11-03 09:18:49Z daan $

ca = zeros(s,1);
cxi1 = zeros(s,1);
cxi2 = zeros(s,1);

rho_xi = hoi_asy_rho_xi(fxi, gxi, s-1);
rho_xi_d = hoi_asy_rho_xi_deriv(fxi, gxi, s-1);
rho_a = hoi_asy_rho(fa, ga, rho_xi, s-1);

for m=1:s
    cxi1(m) = (1i)^(m-1) * rho_xi(m);
    ca(m) = -(1i)^m * (rho_a(m)-rho_xi(m)) / ga(2);
    cxi2(m) = -(1i)^m * rho_xi_d(m) / gxi(3);
end
