% RADAU_JACOBI  Gauss-Radau quadrature rule for Jacobi weight 
% function.
%
%    xw=RADAU_JACOBI(n,iopt,a,b) generates the (n+1)-point Gauss-Radau
%    rule for the Jacobi weight function on [-1,1] with parameters
%    a and b. If iopt=1, the fixed node is at the left endpoint,
%    otherwise at the right endpoint, of the interval [-1,1]. The 
%    n+1 nodes are stored in the first column, the n+1 weights in
%    the second column, of the (n+1)x2 array xw.  The call xw=
%    RADAU_JACOBI(n,iopt,a) is the same as xw=RADAU_JACOBI(n,iopt,
%    a,a), and xw=RADAU_JACOBI(n,iopt) the same as xw=RADAU_JACOBI
%    (n,iopt,0,0).
%
%    REFERENCE: W. Gautschi, ``Gauss-Radau formulae for Jacobi and
%    Laguerre weight functions'', Math. Comput. Simulation 54
%    (2000), 403-412.
%
function xw=radau_jacobi(N,iopt,a,b)
if nargin<2, a=0; end; if nargin<3, b=a; end
ab=r_jacobi(N+1,a,b);
if iopt==1
  ab(N+1,1)=-1+2*N*(N+a)/((2*N+a+b)*(2*N+a+b+1));
else
  ab(N+1,1)=1-2*N*(N+b)/((2*N+a+b)*(2*N+a+b+1));
end
xw=gauss(N+1,ab);
