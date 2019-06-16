% CHRI3  Modification by a simplified symmetric quadratic factor.
%
%    Given a symmetric weight function w(t) through the first n+2
%    recurrence coefficients ab0 of its orthogonal polynomials,
%    ab=CHRI3(n,ab0,y) generates the first n recurrence
%    coefficients of the orthogonal polynomials relative to the
%    modified weight function (t^2+y^2)w(t). The alpha- and beta-
%    coefficients of the original weight function are to be 
%    provided in the first and second column of the (n+1)x2 array
%    ab0; those of the modified weight function are returned in the
%    first and second column of the nx2 array ab.
%
function ab=chri3(N,ab0,y)
N0=size(ab0,1); if N0<N+2, error('input array ab0 too short'), end
r=zeros(N+2,1);
r(1)=y; r(2)=y+ab0(2,2)/y; r(3)=y+ab0(3,2)/r(2);
ab(1,1)=0; ab(1,2)=ab0(1,2)*(ab0(2,2)+y^2);
if N==1, return, end
for k=2:N
  r(k+2)=y+ab0(k+2,2)/r(k+1);
  ab(k,1)=0; ab(k,2)=ab0(k,2)*r(k+1)/r(k-1);
end
