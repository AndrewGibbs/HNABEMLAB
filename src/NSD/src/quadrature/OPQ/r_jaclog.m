% R_JACLOG  Recurrence coefficients for an algebraic/logarithmic
% weight function.
%
%    ab=R_JACLOG(n,a) generates the first n recurrence coefficients
%    for monic polynomials orthogonal with respect to the weight
%    function w(t)=t^a*log(1/t) on [0,1]. The call ab=R_JACLOG(n) 
%    is the same as ab=R_JACLOG(n,0).
%
function ab=r_jaclog(N,a)
if nargin<2, a=0; end;  
if((N<=0)|(a<=-1)) error('parameter(s) out of range'), end
abj=r_jacobi01(2*N);
mom=mm_log(2*N,a);
ab=chebyshev(N,mom,abj);
