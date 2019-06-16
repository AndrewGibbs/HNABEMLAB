% RADAU_LAGUERRE  Gauss-Radau quadrature rule for Laguerre weight 
% function.
%
%    xw=RADAU_LAGUERRE(n,a) generates the (n+1)-point Gauss-Radau
%    rule for the Laguerre weight function on [0,Inf] with parameter a.
%    The n+1 nodes are stored in the first column, the n+1 weights
%    in the second column, of the (n+1)x2 array xw.  The call xw=
%    RADAU_LAGUERRE(n) is the same as xw=RADAU_LAGUERRE(n,0).
%
%    REFERENCE: W. Gautschi, ``Gauss-Radau formulae for Jacobi and
%    Laguerre weight functions'', Math. Comput. Simulation 54
%    (2000), 403-412.
%
function xw=radau_laguerre(N,a)
if nargin<2, a=0; end
ab=r_laguerre(N+1,a);
ab(N+1,1)=N;
xw=gauss(N+1,ab);
