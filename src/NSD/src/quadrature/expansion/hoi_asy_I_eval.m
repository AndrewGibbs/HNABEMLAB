function z = hoi_asy_I_eval(s, omega, ca, cb, ga, gb)
%function z = hoi_asy_I_eval(s, omega, ca, cb, ga, gb)
%
%   Evaluate the asymptotic expansion returned by hoi_asy_I.

% $Id: hoi_asy_I_eval.m 2 2006-11-03 09:18:49Z daan $

z = 0;
O = 1;
for m=1:s
    c = exp(1i*omega*gb(1))*cb(m) - exp(1i*omega*ga(1))*ca(m);
    O = O*omega;
    z = z + c/O;
end
