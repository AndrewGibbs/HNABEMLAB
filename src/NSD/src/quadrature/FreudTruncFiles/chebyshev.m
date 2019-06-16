% CHEBYSHEV Modified Chebyshev algorithm.
%
%    Given a weight function w encoded by its first 2n modified
%    moments, stored in the (row) vector mom, relative to monic 
%    polynomials defined by the (2n-1)x2 array abm of their
%    recurrence coefficients, [ab,normsq]=CHEBYSHEV(n,mom,abm)
%    generates the array ab of the first n recurrence coefficients
%    of the orthogonal polynomials for the weight function w, and 
%    the vector normsq of their squared norms. The n alpha-
%    coefficients are stored in the first column, the n beta-
%    coefficients in the second column, of the nx2 array ab. The
%    call [ab,normsq]=CHEBYSHEV(n,mom) does the same, but using the 
%    classical Chebyshev algorithm. If n is larger than the sizes
%    of mom and abm warrant, then n is reduced accordingly.
%
function [ab,normsq]=chebyshev(N,mom,abm)
if N<=0, error('N out of range'), end
if N>size(mom,2)/2, N=size(mom,2)/2; end
if nargin<3, abm=zeros(2*N-1,2); end
if N>(size(abm,1)+1)/2, N=(size(abm,1)+1)/2; end
ab(1,1)=abm(1,1)+mom(2)/mom(1); ab(1,2)=mom(1);
if N==1, normsq(1)=mom(1); return, end
sig(1,1:2*N)=0; sig(2,:)=mom(1:2*N);
for n=3:N+1
  for m=n-1:2*N-n+2
    sig(n,m)=sig(n-1,m+1)-(ab(n-2,1)-abm(m,1))*sig(n-1,m) ...
      -ab(n-2,2)*sig(n-2,m)+abm(m,2)*sig(n-1,m-1);
  end
  ab(n-1,1)=abm(n-1,1)+sig(n,n)/sig(n,n-1)-sig(n-1,n-1)/ ...
    sig(n-1,n-2);
  ab(n-1,2)=sig(n,n-1)/sig(n-1,n-2);
end
for n=1:N, normsq(n)=sig(n+1,n); end; normsq=normsq'; 
