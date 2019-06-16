function [x,w] = sinc_quadrature(N,h)
%function [x,w] = sinc_quadrature(N,h)
%
%   Construct a sinc quadrature rule with 2N+1 points.

j = (-N:N)';
x=(exp(j*h)-1)./(exp(j*h)+1);
w=2*h.*exp(j*h)./ (1+exp(j*h)).^2;

