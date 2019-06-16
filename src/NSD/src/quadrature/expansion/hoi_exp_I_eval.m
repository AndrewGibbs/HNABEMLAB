function z = hoi_exp_I_eval(s, omega, ca, cb, ga, gb)
%function z = hoi_exp_I_eval(s, omega, ca, cb, ga, gb)
%
%   Evaluate the asymptotic expansion returned by hoi_exp_I.

% $Id: hoi_exp_I_eval.m 2 2006-11-03 09:18:49Z daan $

z = hoi_exp_F_eval(s, omega, cb, gb) - hoi_exp_F_eval(s, omega, ca, ga);
