function z = hoi_exp_F_spr_eval(s, r, omega, ca, ga)
%function z = hoi_exp_F_spr_eval(s, r, omega, ca, ga)
%
%   Evaluate the asymptotic expansion returned by hoi_expf_F_spr.

% $Id$

z = 0;
O = 1;
for m=1:s
    c = exp(1i*omega*ga(1))*ca(m);
    O = O*omega^(1/(r+1));
    z = z + c/O;
end
