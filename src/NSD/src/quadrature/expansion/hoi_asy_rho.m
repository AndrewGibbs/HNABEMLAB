function rho = hoi_asy_rho(fdata, gdata, rho_xi, s)
%function rho = hoi_asy_rho(fdata, gdata, rho_xi, s)
%
%   Compute the functions rho for the asymptotic expansion at a regular
%   point.
%   These require the derivatives of f up to order s, and derivatives of g
%   up to order s+1.

% $Id: hoi_asy_rho.m 2 2006-11-03 09:18:49Z daan $

rho = zeros(s+1, 1);

a0 = fdata(1);
b1 = gdata(2);

rho(1) = a0;

if s >= 1
    a1 = fdata(2);
    b2 = gdata(3);
    rxi0 = rho_xi(1);
    rho(2) = a1/b1-(a0-rxi0)*b2/b1^2;
end

if s >= 2
    a2 = fdata(3);
    b3 = gdata(4);
    rxi1 = rho_xi(2);
    rho(3) = (a2*b1^2-3*a1*b2*b1+3*b2^2*a0-3*b2^2*rxi0-b3*b1*a0+b3*b1*rxi0+b2*rxi1*b1^2)/b1^4;
end

if s >= 3
    a3 = fdata(4);
    b4 = gdata(5);
    rxi2 = rho_xi(3);
    rho(4) = (a3*b1^3-6*a2*b2*b1^2+15*a1*b2^2*b1-4*a1*b3*b1^2-15*b2^3*a0+15*b2^3*rxi0+10*b2*b3*b1*a0-10*b2*b3*b1*rxi0-b4*b1^2*a0+b4*b1^2*rxi0-3*b2^2*rxi1*b1^2+b3*b1^3*rxi1+b2*rxi2*b1^4)/b1^6;
end

if s >= 4
    a4 = fdata(5);
    b5 = gdata(6);
    rxi3 = rho_xi(4);
    rho(5) = (10*b3^2*b1^2*a0-b5*b1^3*a0-10*b3^2*b1^2*rxi0+b4*b1^4*rxi1+15*b2^3*rxi1*b1^2+105*b2^4*a0-105*b2^4*rxi0-3*b2^2*rxi2*b1^4+b5*b1^3*rxi0+a4*b1^4+b3*b1^5*rxi2+60*a1*b2*b3*b1^2+15*b2*b4*b1^2*a0-15*b2*b4*b1^2*rxi0-10*a3*b2*b1^3+45*a2*b2^2*b1^2-10*a2*b3*b1^3-105*a1*b2^3*b1-5*a1*b4*b1^3-105*b2^2*b3*b1*a0+105*b2^2*b3*b1*rxi0+b2*rxi3*b1^6-10*b2*b3*b1^3*rxi1)/b1^8;
end

if s >= 5
    warningdata('s too high in hoi_asy_rho');
end

