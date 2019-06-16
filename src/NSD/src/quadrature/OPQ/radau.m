% RADAU Gauss-Radau quadrature rule.
%
%    Given a weight function w encoded by the (n+1)x2 array ab of
%    the first n+1 recurrence coefficients for the associated
%    orthogonal polynomials, the first column of ab containing the
%    n+1 alpha-coefficients and the second column the n+1 beta-
%    coefficients, the call xw=RADAU(n,ab,end0) generates the 
%    nodes and weights xw of the (n+1)-point Gauss-Radau 
%    quadrature rule for the weight function w having a prescribed
%    node end0 (typically at one of the end points of the support
%    interval of w, or outside thereof). The n+1 nodes, in 
%    increasing order, are stored in the first column, the n+1 
%    corresponding weights in the second column, of the (n+1)x2 
%    array xw.
%
%    For Jacobi weight functions, see also RADAU_JACOBI.
%
function xw=radau(N,ab,end0)
N0=size(ab,1); if N0<N+1, error('input array ab too short'), end
ab0=ab;
p0=0; p1=1;
for n=1:N
  pm1=p0; p0=p1;
  p1=(end0-ab0(n,1))*p0-ab0(n,2)*pm1;
end
ab0(N+1,1)=end0-ab0(N+1,2)*p0/p1;
xw=gauss(N+1,ab0);
