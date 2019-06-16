function rho = hoi_asy_rho_xi(fdata, gdata, s)
%function rho = hoi_asy_rho_xi(fdata, gdata, s)
%
%   Compute the functions rho for the asymptotic expansion at a stationary
%   point xi of order 1.
%   These require the derivatives of f up to order 2s, and derivatives of g
%   up to order 2s+1.

% $Id: hoi_asy_rho_xi.m 2 2006-11-03 09:18:49Z daan $

rho = zeros(s+1, 1);

a0 = fdata(1);

rho(1) = a0;

if s >= 1
    a1 = fdata(2);
    a2 = fdata(3);
    c2 = gdata(3);
    c3 = gdata(4);
    rho(2) = (c2*a2-a1*c3)/(2*c2^2);
end

if s >= 2
    a3 = fdata(4);
    a4 = fdata(5);
    c4 = gdata(5);
    c5 = gdata(6);
    rho(3) = -1/24*(10*c2^2*a3*c3-15*a2*c2*c3^2+15*a1*c3^3-16*a1*c3*c4*c2+6*c2^2*a2*c4+3*c2^2*a1*c5-3*c2^3*a4)/c2^5;
end

if s >= 3
    a5 = fdata(6);
    a6 = fdata(7);
    c6 = gdata(7);
    c7 = gdata(8);
    rho(4) = -1/48*(-c2^5*a6+60*c2^2*a3*c3^3+3*c2^4*a2*c6+90*a1*c3^5-10*c2^3*a1*c6*c3+125*c2^2*a2*c4*c3^2-25*c2^3*a4*c3^2+7*c2^4*a5*c3-15*c2^3*a1*c4*c5+66*c2^2*a1*c4^2*c3-28*c2^3*a2*c3*c5+53*c2^2*a1*c3^2*c5-185*c2*a1*c4*c3^3-50*c2^3*a3*c3*c4+8*c2^4*a4*c4+c2^4*a1*c7-16*c2^3*a2*c4^2-90*c2*a2*c3^4+7*c2^4*a3*c5)/c2^8;
end

if s >= 4
    a7 = fdata(8);
    a8 = fdata(9);
    c8 = gdata(9);
    c9 = gdata(10);
    rho(5) = 1/5760*(-49725*a1*c3^7-33150*c2^2*a3*c3^5-48585*c2^2*a1*c5*c3^4-121400*c2^2*a1*c4^2*c3^3+33810*c2^3*a2*c5*c3^3+62300*c2^3*a2*c4^2*c3^2-17800*c2^4*a3*c3*c4^2-14490*c2^4*a3*c5*c3^2-2030*c2^4*a1*c3^2*c7-6780*c2^4*a2*c3^2*c6+151350*c2*a1*c3^5*c4+59100*c2^3*a3*c3^3*c4-17800*c2^4*a4*c4*c3^2+1500*c2^5*a2*c4*c6+900*c2^5*a2*c3*c7-4305*c2^4*a1*c3*c5^2+678*c2^5*a1*c5*c6-5340*c2^4*a1*c4^2*c5+11610*c2^3*a1*c3^3*c6+240*c2^5*a1*c8*c3+3360*c2^5*a5*c3*c4+480*c2^5*a1*c4*c7+3360*c2^5*a4*c3*c5+21760*c2^3*a1*c4^3*c3+3360*c2^5*a3*c4*c5+2260*c2^5*a3*c3*c6-118200*c2^2*a2*c4*c3^4+1130*c2^5*a6*c3^2-15*c2^6*a1*c9-60*c2^6*a2*c8-20160*c2^4*a2*c3*c4*c5-7120*c2^4*a1*c4*c6*c3+52450*c2^3*a1*c3^2*c4*c5-180*c2^6*a3*c7+14775*c2^3*a4*c3^4-180*c2^6*a7*c3-4830*c2^4*a5*c3^3-300*c2^6*a6*c4-378*c2^6*a5*c5-300*c2^6*a4*c6+49725*c2*a2*c3^6+945*c2^5*a2*c5^2+1980*c2^5*a4*c4^2-3960*c2^4*a2*c4^3+15*c2^7*a8)/c2^11;
end

% if s >= 5
%     a9 = fdata(10);
%     a10 = fdata(11);
%     c10 = gdata(11);
%     c11 = gdata(12);
%     rho(6) = 0;
% end

if s >= 5
    warningdata('s too high in hoi_asy_rho_xi');
end

