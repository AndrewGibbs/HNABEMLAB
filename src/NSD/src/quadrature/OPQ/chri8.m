% CHRI8 Modification by a quadratic divisor
%
%    This routine generates the first N recurrence coefficients of 
%    the orthogonal polynomials relative to the modified weight function 
%    w(t)/(t^2-x^2), where x is real outside the support of w and
%    the weight function w is given by the recurrence coefficients
%    of its orthogonal polynomials. If iopt=1, it computes the desired
%    recurrence coefficients by means of ratios of appropriate Cauchy
%    integrals. These are computed to a relative accuracy eps0 by the
%    continued fraction algorithm implemented in the routine cauchy.m.
%    This calls for the input parameters nu0, numax, and the input
%    array ab0 of size (numax)x2 containing the first numax recurrence
%    coefficients of the orthogonal polynomials relative to the weight
%    function w. The desired coefficients are returned in the Nx2
%    output array ab. The output parameter nu is inherited from the
%    routine cauchy.m. If iopt~=1, the routine uses forward recursion
%    to compute the required ratios of Cauchy integrals. In that case,
%    the routine needs the values of rho0p=rho0(x) and rho0m=rho0(-x)
%    of the first Cauchy integral, and it returns nu=0.
%
function [ab,nu]=chri8(N,ab0,x,eps0,nu0,numax,rho0p,rho0m,iopt)
if N<1, error('N out of range'), end
N0=size(ab0,1); rp=zeros(N+1,1); rm=zeros(N+1);
s=zeros(N,1); 
if N>=2, t=zeros(N-1,1); end
if iopt==1
  if N0<numax, error('input array ab0 too short'), end
  [rhop,rp,nu]=cauchy(N,ab0,x,eps0,nu0,numax);
  [rhom,rm,num]=cauchy(N,ab0,-x,eps0,nu0,numax);
  rho0p=rhop(1); rho0m=rhom(1);
else
  if N0<N, error('input array too short'), end
  nu=0;
  rp(1)=rho0p; rm(1)=rho0m;
  for k=1:N
    rp(k+1)=x-ab0(k,1)-ab0(k,2)/rp(k);
    rm(k+1)=-x-ab0(k,1)-ab0(k,2)/rm(k);
  end
end
for k=1:N
  s(k)=-(rm(k+1)/rp(k)-rp(k+1)/rm(k))/(1/rp(k)-1/rm(k));
  if (k>=2)
    t(k)=(rm(k+1)-rp(k+1))/(1/rp(k)-1/rm(k));
  end
end
ab(1,1)=x*(rho0p+rho0m)/(rho0p-rho0m);
ab(1,2)=-(rho0p-rho0m)/(2*x);
if N==1, break, end
ab(2,1)=ab0(2,1)-s(2)+s(1);
ab(2,2)=ab0(2,2)+s(1)*(ab0(1,1)-ab(2,1))-t(2);
if N==2, break, end
ab(3,1)=ab0(3,1)-s(3)+s(2);
ab(3,2)=ab0(3,2)+s(2)*(ab0(2,1)-ab(3,1))-t(3)+t(2);
if N==3, break, end
for k=4:N
  ab(k,1)=ab0(k,1)-s(k)+s(k-1);
  ab(k,2)=ab0(k-2,2)*t(k-1)/t(k-2);
end
