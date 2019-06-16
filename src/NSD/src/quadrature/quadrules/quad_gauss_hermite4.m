function [x,w] = quad_gauss_hermite4(N)
%function [x,w] = quad_gauss_hermite4(N)
%
%   Return the quadrature rule with N points for the weight function
%   exp(-x^4) on (-infty,infty).
%
%   References:
%   W. Gautschi, ``Orthogonal Polynomials: Computation and Approximation'',
%       Clarendon Press, Oxford, 2004.

% $Id$

beta = [1.8128049541109541560, .33798912003364236450, .40167965976351735858, .50510423234482229782, .57805815033171132110, .64676738204724497038, .70786315090515246154, .76442312605207732341, .81702175201098194520, .86647036419002229799, .91324989944000747280, .95775608488417898910, 1.0002887465597799483, 1.0410891789282270540, 1.0803527252385584741, 1.1182407644978825284, 1.1548882639750164576]';
alpha = zeros(size(beta));

[x,w] = compute_gauss(N, alpha, beta);

