% CHRI2  Modification by a quadratic factor.
%
%    Given a weight function w(t) through the first n+2 recurrence
%    coefficients ab0 of its orthogonal polynomials, ab=CHRI2(n,
%    ab0,x,y) generates the first n recurrence coefficients 
%    of the orthogonal polynomials relative to the modified weight
%    function ((t-x)^2+y^2)w(t). The alpha- and beta-coefficients
%    of the original weight function are to be provided in the 
%    first and second column of the (n+2)x2 array ab0; those of the
%    modified weight function are returned in the first and second
%    column of the nx2 array ab.
%
function ab=chri2(N,ab0,x,y)
N0=size(ab0,1); if N0<N+2, error('input array ab0 too short'), end
r=zeros(N+2,1);
z=complex(x,y);
r(1)=z-ab0(1,1); r(2)=z-ab0(2,1)-ab0(2,2)/r(1);
r(3)=z-ab0(3,1)-ab0(3,2)/r(2);
ab(1,1)=ab0(3,1)+real(r(3))+(imag(r(3))/imag(r(2)))*real(r(2))...
  -(real(r(2))+(imag(r(2))/imag(r(1)))*real(r(1)));
ab(1,2)=ab0(1,2)*(ab0(2,2)+abs(r(1))^2);
if N==1, return, end
for k=1:N-1
  r(k+3)=z-ab0(k+3,1)-ab0(k+3,2)/r(k+2);
  ab(k+1,1)=ab0(k+3,1)+real(r(k+3))+(imag(r(k+3))/imag(r(k+2)))...
    *real(r(k+2))-(real(r(k+2))+(imag(r(k+2))/imag(r(k+1)))*real(r(k+1)));
  ab(k+1,2)=ab0(k+1,2)*(imag(r(k+2))*imag(r(k))/imag(r(k+1))^2)...
    *abs(r(k+1)/r(k))^2;
end
