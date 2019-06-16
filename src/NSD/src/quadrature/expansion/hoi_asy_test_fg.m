function [f, g, gd,  fdata, gdata] = hoi_asy_test_fg(x)
%function [f, g, gd, fdata, gdata] = hoi_asy_test_fg(x)
%
%   Return test functions f and g and their derivatives in the vectors
%   fdata and gdata. There are no stationary points on [0,1].

% $Id: hoi_asy_test_fg.m 2 2006-11-03 09:18:49Z daan $


f = @test_f;
g = @test_g;
gd = @test_gd;
fdata = [x.^5 + x.^4 + x.^3 + x.^2 + x + 1 ...
    5*x.^4 + 4*x.^3+3*x.^2+2*x+1 ...
    20*x.^3 + 12*x.^2 + 6*x+2 ...
    60*x.^2 + 24*x + 6 ...
    120*x + 24];
gdata = [x.^6 + x.^5 + x.^4 + x.^3 + x.^2 + x + 1 ...
    6*x.^5 + 5*x.^4 + 4*x.^3+3*x.^2+2*x+1 ...
    30*x.^4 + 20*x.^3 + 12*x.^2 + 6*x+2 ...
    120*x.^3 + 60*x.^2 + 24*x + 6 ...
    360*x.^2 + 120*x + 24 ...
    720*x + 120];


function z = test_f(x)
z = x.^5 + x.^4 + x.^3 + x.^2 + x + 1;


function z = test_g(x)
z = x.^6 + x.^5 + x.^4 + x.^3 + x.^2 + x + 1;

function z = test_gd(x)
z = 6*x.^5 + 5*x.^4 + 4*x.^3 + 3*x.^2 + 2*x + 1;

