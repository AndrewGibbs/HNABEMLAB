function z = hoi_asy_I_sp1b_eval(s, omega, ca, cxi1, cxi2, ga, gxi, mu0)
%function z = hoi_asy_I_sp1b_eval(s, omega, ca, cxi1, cxi2, ga, gxi, mu0)
%
%   Evaluate the asymptotic expansion returned by hoi_asy_I_sp1b.

% $Id: hoi_asy_I_sp1b_eval.m 2 2006-11-03 09:18:49Z daan $

z = 0;
O = 1;
for m=1:s
    c1 = mu0*cxi1(m)/O;
    O = O*omega;
    c2 = exp(1i*omega*gxi(1))*cxi2(m)/O - exp(1i*omega*ga(1))*ca(m)/O;
    z = z + c1 + c2;
end
