% ACONDG Absolute condition number of G_n.
%
%    This computes the absolute condition number of the map G_n
%    from the first 2n normalized modified moments of a measure
%    dlambda to the n-point Gauss quadrature formula. The measure
%    is specified via the first n recursion coefficients, stored
%    in the nx2 array ab, of the polynomials orthogonal with respect
%    to the measure.  The orthogonal polynomials defining the
%    modified moments are specified via the (2n)x2 array ab0 of
%    their first 2n recursion coefficients.
%
function ac=acondG(n,ab,ab0)
g=zeros(2*n,1);
xw=gauss(n,ab); xw0=gauss(2*n,ab0); 
for nu=1:2*n
  g(nu)=g_n(n,xw0(nu,1),xw);
end
ac=sqrt(sum(xw0(:,2).*g));
