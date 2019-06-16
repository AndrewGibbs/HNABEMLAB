% R_ELLIPTIC Recurrence coefficients for an elliptic weight function.
%
%    ab=R_ELLIPTIC(n,k2) generates the first n recurrence coefficients
%    for monic polynomials orthogonal with respect to the weight
%    function w(t)=[(1-k2*t^2)*(1-t^2)]^{-1/2} on [-1,1], where
%    0<k2<1.
%
function ab=r_elliptic(N,k2)
if((N<=0)|(k2<0)|(k2>=1)), error('parameter(s) out of range'), end
abm=r_jacobi(2*N-1,-1/2);
mom=mm_ell(N,k2);
ab=chebyshev(N,mom,abm);
