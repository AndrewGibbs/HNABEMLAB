% R_MEIXNER_POLLACZEK Recurrence coefficients for monic 
% Meixner-Pollaczek polynomials.
%
%    ab=R_MEIXNER_POLLACZEK(n,lambda,phi) generates the first n 
%    recurrence coefficients for the Meixner-Pollaczek polynomials
%    with parameters lambda and phi. These are orthogonal on 
%    [-Inf,Inf] relative to the weight function w(t)=(2 pi)^(-1)
%    exp((2 phi-pi)t) |gamma(lambda+i t)|^2. The n alpha-
%    coefficients are stored in the first column, the n beta-
%    coefficients in the second column, of the nx2 array ab. The
%    call ab=R_MEIXNER_POLLACZEK(n,lambda) is the same as ab=
%    R_MEIXNER_POLLACZEK(n,lambda,pi/2).
%
function ab=r_meixner_pollaczek(N,lambda,phi)
if nargin<3, phi=pi/2; end
if ((N<=0)|(lambda<=0)|(phi<=0)|(phi>=pi)) 
  error('parameter(s) out of range') 
end
n=1:N; sinphi=sin(phi); lam2=2*lambda;
if sinphi==1
  ab(:,1)=0;
else
  a(:,1)=-(n+lambda-1)/tan(phi);
end
ab(:,2)=(n-1).*(n+lam2-2)/(4*sinphi^2);
ab(1,2)=gamma(lam2)/(2*sinphi)^lam2;
