% GCONDG Global condition number of G_n.
%
%    This computes the global condition number of the map G_n from
%    the first 2n ordinary moments of a measure dlambda to the n-point 
%    Gauss quadrature formula. The measure is specified via the first 
%    n recursion coefficients, stored in the nx2 array ab, of the
%    polynomials orthogonal with respect to the measure, whereas 
%    the moments are input through the (2n)-vector mom.
%
function gc=gcondG(n,ab,mom)
xw=gauss(n,ab);
Linv=diag([ones(n,1);1./xw(:,2)]);
gam=[xw(:,2);xw(:,1)];
gc=norm(mom,inf)*norm(Linv*Tinv(n,xw),inf)/norm(gam,inf);
