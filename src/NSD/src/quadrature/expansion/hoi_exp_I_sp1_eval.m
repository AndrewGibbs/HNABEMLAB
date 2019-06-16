function z = hoi_exp_I_sp1_eval(s, n, omega, ca, cb, cxi, ga, gb, gxi)
%function z = hoi_exp_I_sp1_eval(s, n, omega, ca, cb, cxi, ga, gb, gxi)
%
%   Evaluate the asymptotic expansion returned by hoi_exp_I_sp1.

% $Id: hoi_exp_I_sp1_eval.m 2 2006-11-03 09:18:49Z daan $

z = hoi_exp_F_eval(s, omega, cb, gb) ...
    + hoi_exp_F_sp1_eval(n, omega, cxi, gxi) ...
    - hoi_exp_F_eval(s, omega, ca, ga);
