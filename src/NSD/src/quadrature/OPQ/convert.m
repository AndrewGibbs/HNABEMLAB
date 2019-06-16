% CONVERT Conversion algorithm
%
%    Given the N+1 coefficients, input in the (N+1)x1 array c0, of the
%    expansion of a polynomial p of degree N in orthogonal polynomials
%    defined by the (N+1)x2 array ab0 of their recurrence coefficients,
%    the call c=CONVERT(N,c0,ab0,ab) generates and stores in the
%    (N+1)x1 array c, the N+1 coefficients of the same polynomial p
%    expanded in orthogonal polynomials defined by the (N+1)x2 array ab
%    of their recurrence coefficients. The choice ab0=0 is permitted,
%    in which case the given expansion is a power series expansion.
%
function c=convert(N,c0,ab0,ab)
sig=zeros(N+3);
sig(2,2)=ab(1,2); 
for m=1:N
  for n=1:m+1
    sig(n+1,m+2)=sig(n+2,m+1)+(ab(n,1)-ab0(m,1))*sig(n+1,m+1)...
      +ab(n,2)*sig(n,m+1)-ab0(m,2)*sig(n+1,m);
  end
end
for k=1:N+1
  c(k)=sum(sig(k+1,k+1:N+2).*c0(k:N+1))/sig(k+1,k+1);
end
