% LOBATTO Gauss-Lobatto quadrature rule.
%
%    Given a weight function w encoded by the (n+2)x2 array ab of
%    the first n+2 recurrence coefficients for the associated
%    orthogonal polynomials, the first column of ab containing the
%    n+2 alpha-coefficients and the second column the n+2 beta-
%    coefficients, the call xw=LOBATTO(n,ab,endl,endr) generates 
%    the nodes and weights xw of the (n+2)-point Gauss-Lobatto 
%    quadrature rule for the weight function w having two
%    prescribed nodes endl, endr (typically the left and right end
%    points of the support interval of w, or points to the left
%    resp. to the right therof). The n+2 nodes, in increasing 
%    order, are stored in the first column, the n+2 corresponding
%    weights in the second column, of the (n+2)x2 array xw.
%        
%    For Jacobi weight functions, see also LOBATTO_JACOBI.
%
function xw=lobatto(N,ab,endl,endr)
N0=size(ab,1); if N0<N+2, error('input array ab too short'), end
ab0=ab;
p0l=0; p0r=0; p1l=1; p1r=1;
for n=1:N+1
  pm1l=p0l; p0l=p1l; pm1r=p0r; p0r=p1r;
  p1l=(endl-ab0(n,1))*p0l-ab0(n,2)*pm1l;
  p1r=(endr-ab0(n,1))*p0r-ab0(n,2)*pm1r;
end
det=p1l*p0r-p1r*p0l;
ab0(N+2,1)=(endl*p1l*p0r-endr*p1r*p0l)/det;
ab0(N+2,2)=(endr-endl)*p1l*p1r/det;
xw=gauss(N+2,ab0);
