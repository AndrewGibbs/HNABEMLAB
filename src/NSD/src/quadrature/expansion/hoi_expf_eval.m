function z = hoi_expf_eval(c, g0, omega, r)
%function z = hoi_expf_eval(c, g0, omega, r)
%
%   Evaluate the given asymptotic expansion.

% $Id$

z = 0;

s = length(c);
for i=1:s
    z = z + c(i) * omega^(-i/(r+1));
end
z = exp(1i*omega*g0)*z;
