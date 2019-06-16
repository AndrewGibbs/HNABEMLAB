function [c1,c2] = hoi_exp_F_sp1(n, xi, fxi, gxi)
%function [c1,c2] = hoi_exp_F_sp1(n, xi, fxi, gxi)
%
%   Compute n terms of the asymptotic expansion of F at a stationary point
%   xi of order 1.
%   These n terms require n-1 derivatives of f if n is even, and n
%   derivatives of f if n is odd. They require n+1 derivatives of g.

% $Id: hoi_exp_F_sp1.m 2 2006-11-03 09:18:49Z daan $


[d1,d2] = hoi_exp_F0_sp1(n, xi, gxi);

s = ceil(n/2);
rho_xi = hoi_asy_rho_xi(fxi, gxi, s-1);
rho_xi_d = hoi_asy_rho_xi_deriv(fxi, gxi, s-1);
cxi1 = zeros(s,1);
cxi2 = zeros(s,1);
for m=1:s
    cxi1(m) = (1i)^(m-1) * rho_xi(m);
    cxi2(m) = -(1i)^m * rho_xi_d(m) / gxi(3);
end

c1 = zeros(n, 1);
c2 = zeros(n, 1);

for i=1:n
    % find the coefficient in omega^(-i/2)
    e1 = 0;
    e2 = 0;
    for j=1:i
        for l=0:floor(i/2)
            if (j+2*l) == i
                e1 = e1 + d1(j)*cxi1(l+1);
                e2 = e2 + d2(j)*cxi1(l+1);
            end
        end
    end
    c1(i) = e1;
    c2(i) = e2;
    if (mod(i,2) == 0)
        c1(i) = c1(i) + cxi2(i/2);
        c2(i) = c2(i) + cxi2(i/2);
    end
end
