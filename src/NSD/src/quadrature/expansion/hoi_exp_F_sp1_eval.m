function z = hoi_exp_F_sp1_eval(n, omega, cxi, gxi)
%function z = hoi_exp_F_sp1_eval(n, omega, cxi, gxi)
%
%   Evaluate the asymptotic expansion returned by hoi_exp_F_sp1.

% $Id: hoi_exp_F_sp1_eval.m 2 2006-11-03 09:18:49Z daan $

z = 0;
O = 1;
for m=1:n
    c = exp(1i*omega*gxi(1))*cxi(m);
    O = O*sqrt(omega);
    z = z + c/O;
end
