function [x,w] = quad_gauss_jacobi(N, alpha, beta)
%function [x,w] = quad_gauss_jacobi(N, alpha, beta)
%
%   Return the Gaussian quadrature with N points in the interval
%   [-1,1] for the weight function w(x)=(1-t)^alpha (1+t)^beta.
%
%   References:
%   W. Gautschi, ``Orthogonal Polynomials: Computation and Approximation'',
%       Clarendon Press, Oxford, 2004.

% $Id$

% Use OPQ routine from Gautschi
ab = r_jacobi(N, alpha, beta);
[x,w] = compute_gauss(N, ab(:,1), ab(:,2));
