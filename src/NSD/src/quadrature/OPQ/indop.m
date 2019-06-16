% INDOP Induced orthogonal polynomials.
%
%    Given a weight function w(t) through the first n+m recurrence
%    coefficients ab0 of its orthogonal polynomials, ab=INDOP(n,m,
%    ab0) generates the first n recurrence coefficients of the
%    m-th induced orthogonal polynomial, i.e., the polynomial
%    orthogonal with respect to the weight function pi_m(t)^2w(t),
%    where p_m is the monic polynomial of degree m orthogonal 
%    relative to the weight function w. The alpha- and beta-
%    coefficients of the original weight function are to be provided 
%    in the first and second column of the (n+m)x2 input array ab0;
%    those of the modified weight function are returned in the first 
%    and second column of the nx2 output array ab.
%
function ab=indop(N,m,ab0)
if m<0, error('m out of range'), end
N0=size(ab0,1); if N0<N+m, error('input array ab0 too short'), end
ab=ab0;
if m==0, return, end
zw=gauss(m,ab0);
for imu=1:m
  mi=N+m-imu;
  for n=1:mi+1
    ab1(n,1)=ab(n,1);
    ab1(n,2)=ab(n,2);
  end
  x=zw(imu,1);
  ab=chri7(mi,ab1,x);
end
