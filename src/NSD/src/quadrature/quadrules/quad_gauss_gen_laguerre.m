function [x,w] = quad_gauss_gen_laguerre(N, a)
%function [x,w] = quad_gauss_gen_laguerre(N, a)
%
%   Return the generalized Gauss-Laguerre quadrature rule with N points
%   and parameter a.
%
%   References:
%   W. Gautschi, ``Orthogonal Polynomials: Computation and Approximation'',
%       Clarendon Press, Oxford, 2004.

% $Id$

% Use OPQ routine from Gautschi
ab = r_laguerre(N, a);
[x,w] = compute_gauss(N, ab(:,1), ab(:,2));
