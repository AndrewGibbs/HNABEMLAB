% TINV Inverse of a confluent Vandermonde matrix.
%
%    The call U=TINV(n,x) produces the inverse of a confluent
%    Vandermonde matrix of order 2n whose first n columns are those
%    of an ordinary 2n-th order Vandermonde matrix in the parameters
%    x_1,x_2,...,x_n, which are stored in the input vector x, and 
%    whose remaining n columns are the derivatives of the previous n 
%    columns. The use of the Matlab routine inv is avoided because
%    of ill-condioning. 
%
%    REFERENCE: W. Gautschi, ``On inverses of Vandermonde and confluent
%    Vandermonde matrices. II'', Numer. Math. 5 (1963), 425-430.
%
function U=Tinv(n,x)
x0=zeros(1,n-1);
for lam=1:n
  for nu=1:n-1
    if nu<lam
      x0(nu)=x(nu);
    else
      x0(nu)=x(nu+1);
    end
  end
  tau=zeros(1,2*n+1);
  tau1=zeros(1,2*n-2);
  tau(2)=1;
  for m=0:n-2
    for mu=1:2*m+2
      tau1(mu)=tau(mu+2)+2*x0(m+1)*tau(mu+1)+...
        (x0(m+1)^2)*tau(mu);
    end
    for mu=1:2*m+2, tau(mu+2)=tau1(mu); end
  end
  d=x(lam)-x0;
  s0=sum(1./d); p0=prod(d)^2;
  s=(1+2*x(lam)*s0)/p0; t=2*s0/p0; sgn=1;
  for mu=1:2*n
    U(lam,mu)=sgn*(s*tau(2*n-mu+1)+t*tau(2*n-mu+2));
    sgn=-sgn;
  end
  s=-x(lam)/p0; t=-1/p0;
  for mu=1:2*n
    U(n+lam,mu)=sgn*(s*tau(2*n-mu+1)+t*tau(2*n-mu+2));
    sgn=-sgn;
  end
end
