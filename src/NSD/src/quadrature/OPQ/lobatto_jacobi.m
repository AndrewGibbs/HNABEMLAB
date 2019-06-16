% LOBATTO_JACOBI  Gauss-Lobatto quadrature rule for Jacobi weight 
% function.
%
%    xw=LOBATTO_JACOBI(n,a,b) generates the (n+2)-point Gauss-
%    Lobatto rule for the Jacobi weight function on [-1,1] with
%    parameters a and b. The n+2 nodes are stored in the first
%    column, the n+2 weights in the second column, of the (n+2)x2
%    array xw. The call xw=LOBATTO_JACOBI(n,a) is the same as xw=
%    LOBATTO_JACOBI(n,a,a) and xw=LOBATTO_JACOBI(n) the same as
%    xw=LOBATTO_JACOBI(n,0,0). 
%
%    REFERENCE: W. Gautschi,``High-order Gauss-Lobatto formulae'',
%    Numer. Algorithms 25 (2000), 213-222.
%
function xw=lobatto_jacobi(N,a,b)
if nargin<2, a=0; end; if nargin<3, b=a; end
ab=r_jacobi(N+2,a,b);
ab(N+2,1)=(a-b)/(2*N+a+b+2);
ab(N+2,2)=4*(N+a+1)*(N+b+1)*(N+a+b+1)/((2*N+a+b+1)* ...
(2*N+a+b+2)^2);
xw=gauss(N+2,ab);
