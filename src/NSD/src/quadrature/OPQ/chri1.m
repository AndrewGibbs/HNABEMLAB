% CHRI1 Modification by a linear factor.
%
%    Given a weight function w(t) through the first n+1 recurrence
%    coefficients ab0 of its orthogonal polynomials, ab=CHRI1(n,
%    ab0,z) generates the first n recurrence coefficients for
%    the orthogonal polynomials relative to the modified weight
%    function (t-z)w(t). The alpha- and beta-coefficients of the
%    original weight function are to be provided in the first and
%    second column of the (n+1)x2 array ab0; those of the modified
%    weight function are returned in the first and second column 
%    of the nx2 array ab.
%
function ab=chri1(N,ab0,z)
N0=size(ab0,1); if N0<N+1, error('input array ab0 too short'), end
r=zeros(N+1,1);
r(1)=z-ab0(1,1); r(2)=z-ab0(2,1)-ab0(2,2)/r(1);
ab(1,1)=ab0(2,1)+r(2)-r(1); ab(1,2)=-r(1)*ab0(1,2);
if N==1, return, end
for k=2:N
  r(k+1)=z-ab0(k+1,1)-ab0(k+1,2)/r(k);
  ab(k,1)=ab0(k+1,1)+r(k+1)-r(k);
  ab(k,2)=ab0(k,2)*r(k)/r(k-1);
end
