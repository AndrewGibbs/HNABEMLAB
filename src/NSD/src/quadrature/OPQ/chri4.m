% CHRI4 Modification by a linear divisor
%
%    This routine generates the first N recurrence coefficients of the 
%    orthogonal polynomials relative to the modified weight function 
%    w(t)/(t-z), where z (possibly complex) is outside the support
%    of w and the weight function w is given by the recurrence
%    coefficients of its orthogonal polynomials. If iopt=1, it computes 
%    the desired recurrence coefficients by means of ratios of 
%    appropriate Cauchy integrals. These are computed to a relative 
%    accuracy eps0 by the continued fraction algorithm implemented 
%    in the routine cauchy.m. This calls for the input parameters 
%    nu0, numax, and the input array ab0 of size (numax)x2 containing 
%    the first numax recurrence coefficients of the orthogonal 
%    polynomials relative to the weight function w. The desired 
%    coefficients are returned in the Nx2 output array ab. The output
%    parameter nu is inherited from the routine cauchy.m. If iopt~=1,
%    the routine uses forward recursion to compute the required
%    ratios of Cauchy integrals. In that case, the routine needs the
%    value rho0 of the first Cauchy integral.
%
function [ab,nu]=chri4(N,ab0,z,eps0,nu0,numax,rho0,iopt)
N0=size(ab0,1); r=zeros(N+1,1);
if iopt==1
  if N0<numax, error('input array ab0 too short'), end
  [rho,r,nu]=cauchy(N,ab0,z,eps0,nu0,numax);
else
  if N0<N, error('input array too short'), end
  nu=0;
  r(1)=rho0;
  for k=1:N
    r(k+1)=z-ab0(k,1)-ab0(k,2)/r(k);
  end
end
ab(1,1)=ab0(1,1)+r(2); ab(1,2)=-r(1);
if N==1, return, end
for k=2:N
  ab(k,1)=ab0(k,1)+r(k+1)-r(k);
  ab(k,2)=ab0(k-1,2)*r(k)/r(k-1);
end
