function z = hoi_asy_I_sp1a_eval(s, omega, cxi1, cxi2, cb, gxi, gb, mu0)
%function z = hoi_asy_I_sp1a_eval(s, omega, cxi1, cxi2, cb, gxi, gb, mu0)
%
%   Evaluate the asymptotic expansion returned by hoi_asy_I_sp1a.

% $Id: hoi_asy_I_sp1a_eval.m 2 2006-11-03 09:18:49Z daan $

z = 0;
O = 1;
for m=1:s
    c1 = mu0*cxi1(m)/O;
    O = O*omega;
    c2 = exp(1i*omega*gb(1))*cb(m)/O - exp(1i*omega*gxi(1))*cxi2(m)/O;
    z = z + c1 + c2;
end
