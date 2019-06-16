function z = hoi_exp_I_sp1a_eval(s, n, omega, cxi, cb, gxi, gb)
%function z = hoi_exp_I_sp1a_eval(s, n, omega, cxi, cb, gxi, gb)
%
%   Evaluate the asymptotic expansion returned by hoi_exp_I_sp1a.

% $Id: hoi_exp_I_sp1a_eval.m 2 2006-11-03 09:18:49Z daan $

z = hoi_exp_F_eval(s, omega, cb, gb) ...
    - hoi_exp_F_sp1_eval(n, omega, cxi, gxi);
