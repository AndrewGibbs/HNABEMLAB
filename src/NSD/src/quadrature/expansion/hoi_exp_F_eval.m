function z = hoi_exp_F_eval(s, omega, ca, ga)
%function z = hoi_exp_F_eval(s, omega, ca, ga)
%
%   Evaluate the asymptotic expansion returned by hoi_exp_F.

% $Id: hoi_exp_F_eval.m 2 2006-11-03 09:18:49Z daan $

z = 0;
O = 1;
for m=1:s
    c = exp(1i*omega*ga(1))*ca(m);
    O = O*omega;
    z = z + c/O;
end
