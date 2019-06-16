% BOSE_EINSTEIN_DIFFPOLE  Bose-Einstein integral with difficult pole.
%
function [xw,Ncap,kount,nu1]=bose_einstein_diffpole(N,MM,eta,theta,k,eps0,Nmax,numax)
global mc mp iq idelta irout AB Z ab0 M
%
% Be sure to declare the variable theta global also in the routines
% R_MOD and QUADRAT and make the changes indicated in the routines 
% QUADRAT and FEX32.
%
mc=1; mp=0; iq=1; idelta=2; irout=1; Ncap=0; kount=0; nu1=0;
AB=[0 Inf];
%
% When theta is "large", increase eps0 to avoid an excessively
% large value of Nmax that would be needed for convergence. Also
% increase Nmax appropriatelyr. For example,
%         eps0=1e6*eps; Nmax=800;
%
ab0=r_laguerre(Nmax,k-1);
if MM==0, xw=gauss(N,ab0); return, end
M=MM-1;
if M>1
  for m=1:2:M-1
    Z(m,1)=-1/((m+1)*pi*i); Z(m,2)=1;
    Z(m+1,1)=1/((m+1)*pi*i); Z(m+1,2)=1;
  end
end
[ab1,Ncap,kount]=r_mod(N,ab0);
x1=eta; xa=abs(eta);
t=1; s0=0; s=-1/(k-1); j=0;
while s~=s0
  s0=s; j=j+1;
  if j>100, disp('series for incomplete gamma does not vonverge'), end
  t=-xa*t/j; s=s+t/(j+1-k);
end
hr=-exp(xa)*((pi/sin(pi*k))*(xa^(k-1))-gamma(k)*s);
if M>1
  s=0;
  for n=1:M/2
    p=1;
    for m=1:M/2
      if m~=n
        p=(m^2/(m^2-n^2))*p;
      end
    end
    nu0=10; z=2*n*pi*i;
    [rho,r,nu]=cauchy(0,ab0,z,eps0,nu0,numax);
    if n==1, nu1=nu; end
    s=s+2*n*pi*p*(2*n*pi*hr-real((2*n*pi-x1*i)*rho(1)))/...
      (x1^2+4*n^2*(pi^2));
  end
  hr=s;
end
abmod=chri4(N,ab1,x1,eps0,0,0,hr,0);
abmod(1,2)=-eta*abmod(1,2);
M=M+1;
Z(M,1)=-1/eta; Z(M,2)=1;
xw=gauss_rational(N,abmod);
