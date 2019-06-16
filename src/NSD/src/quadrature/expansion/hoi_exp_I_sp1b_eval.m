function z = hoi_exp_I_sp1b_eval(s, n, omega, ca, cxi, ga, gxi)
%function z = hoi_exp_I_sp1b_eval(s, n, omega, ca, cxi, ga, gxi)
%
%   Evaluate the asymptotic expansion returned by hoi_exp_I_sp1b.

% $Id: hoi_exp_I_sp1b_eval.m 2 2006-11-03 09:18:49Z daan $

z = hoi_exp_F_sp1_eval(n, omega, cxi, gxi) ...
    - hoi_exp_F_eval(s, omega, ca, ga);
