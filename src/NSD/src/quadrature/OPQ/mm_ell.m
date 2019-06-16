% MM_ELL Modified moments for an elliptic weight functions.
%
%    This computes the first 2N Chebyshev moments of the ``elliptic'' weight
%    function w(t)=((1-k^2 t^2)(1-t^2))^{-1/2}. The input parameter k2 has
%    the meaning of k^2.
%
%    REFERENCE:  W. Gautschi, ``On generating orthogonal polynomials'', SIAM J.
%    Sci. Stat. Comput. 3 (1982), 289-317, Example 3.3.
%
function mm=mm_ell(N,k2,eps0)
if nargin<3, eps0=eps; end
if((k2<=0)|(k2>=1)) error('parameter k2 out of range in mm_ell'), end
q=k2/(2-k2+2*sqrt(1-k2));
q1=(1+q^2)/q;
nu=2*N-10; f=ones(1,N); f0=zeros(1,N);
while any(abs(f-f0)>eps0*abs(f))
  nu=nu+10;
  if nu>500, error('nu exceeds 500 in mm_ell'), end
  f0=f; r=0; s=0;
  for n=nu:-1:1
    r=-(n-1/2)/(n*q1+(n+1/2)*r);
    s=r*(2+s);
    if n<=N, rr(n)=r; end
  end
  c0=1/(1+s); f(1)=rr(1)*c0;
  if N>1
    for n=2:N
      f(n)=rr(n)*f(n-1);
    end
  end
end
%
% If the value of nu that yields convergence is of interest,
% it may be printed at this point by adding the statement
% fprintf('nu in mm_ell = %3.0f\n',nu)
%
mm(1)=pi*c0; mm(2)=0; if N==1, return, end
c=2*pi;
for n=3:2:2*N-1
  n1=(n-1)/2;
  c=-c/4;
  mm(n)=c*f(n1); mm(n+1)=0;
end
