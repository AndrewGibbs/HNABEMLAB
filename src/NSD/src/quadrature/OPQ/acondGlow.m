% ACONDGLOW Lower bound for the absolute condition number of G_n.
%
%    This computes a lower bound for the absolute condition number of
%    the map G_n from the first 2n ordinary moments of a measure
%    dlambda to the n-point Gauss quadrature formula. The measure
%    is specified via the first n recursion coefficients, stored 
%    stored in the nx2 array ab, of the polynomials orthogonal with
%    respect to the measure.
%  
function acl=acondGlow(n,ab)
xw=gauss(n,ab);
m=min([1 1/max(xw(:,2))]);
xw=gauss(n,ab);
pnum=prod((1+xw(:,1)).^2);
for nu=1:n
  d=xw(nu,1)-xw(:,1);
  pden(nu)=(1+xw(nu,1))*prod(d(find(d)).^2);
end
acl=m*pnum/min(pden);
