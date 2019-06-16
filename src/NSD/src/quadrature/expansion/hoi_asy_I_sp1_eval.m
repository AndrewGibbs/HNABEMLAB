function z = hoi_asy_I_sp1_eval(s, omega, ca, cb, cxi, ga, gb, mu0)
%function z = hoi_asy_I_sp1_eval(s, omega, ca, cb, cxi, ga, gb, mu0)
%
%   Evaluate the asymptotic expansion returned by hoi_asy_I_sp1.

% $Id: hoi_asy_I_sp1_eval.m 2 2006-11-03 09:18:49Z daan $

z = 0;
O = 1;
for m=1:s
    c1 = mu0*cxi(m)/O;
    O = O*omega;
    c2 = exp(1i*omega*gb(1))*cb(m)/O - exp(1i*omega*ga(1))*ca(m)/O;
    z = z + c1 + c2;
end
