% SPHT Stieltjes-Perron inversion and Hilbert transform.
%
%    If iopt equals 1, this routine tries to compute to an accuracy
%    eps0 or better the weight function w(x) at a point x from the
%    first numax recurrence coefficients of the corresponding
%    orthogonal polynomials to be provided in the numax x 2 array ab.
%    In an attempt to implement the Stieltjes-Perron inversion
%    formula, it generates the Cauchy integral of w for a sequence
%    of n complex arguments approaching x from above and applies
%    the epsilon algorithm to their negative real parts (divided
%    by pi) to accelerate convergence. The accelerated values are
%    returned in the even-numbered columns of the n x (n+1) array E.
%    If iopt is different from 1, the routine tries to compute the
%    Hilbert transform of w at the point x using the same procedure
%    as above, but applying the epsilon algorithm to the negative
%    real parts of the Cauchy integrals.
%    The input parameter nu0 is an estimate of the backward
%    recurrence index to be used in the routine CAUCHY for
%    computing the Cauchy integrals.
%
function [E,w,nu0,nu]=SP(n,x,iopt,ab,nu0,numax,eps0)
N=0; 
w=zeros(n,1); E=zeros(n,n+1);
y=1;
k=0;
while y>0&k<n
  k=k+1;
  y=y/2;
  z=x+i*y; 
  [rho,r,nu]=cauchy(N,ab,z,eps0,nu0,numax);
  if iopt==1
    w(k)=-imag(rho)/pi;
  else
    w(k)=-real(rho);
  end
end
  E=epsalg(n,w);
