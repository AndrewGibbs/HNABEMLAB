% MM_LOG Modified moments for a logarithmic weight function.
%
%    The call mm=MM_LOG(n,a) computes the first n modified moments of the
%    logarithmic weight function w(t)=t^a log(1/t) on [0,1] relative to 
%    shifted Legendre polynomials. 
%
%    REFERENCE:  W. Gautschi,``On the preceding paper `A Legendre polynomial
%    integral' by James L. Blue'', Math. Comp. 33 (1979), 742-743.
%
function mm=mm_log(N,a)
if a<=-1, error('a out of range in mm_log'), end
c=1;
for n=1:N
  if(a-floor(a)==0 & a<n-1)
    p=n-a-1:n+a;
    mm(n)=(-1)^(n-1-a)/prod(p);
    mm(n)=(gamma(a+1))^2*mm(n);
  else
    if n==1
      mm(1)=1/(a+1)^2;
    else 
      k=1:n-1;
      s=1./(a+1+k)-1./(a+1-k);
      p=(a+1-k)./(a+1+k);
      mm(n)=(1/(a+1)+sum(s))*prod(p)/(a+1);
    end
  end
  mm(n)=c*mm(n); c=.5*n*c/(2*n-1);
end
