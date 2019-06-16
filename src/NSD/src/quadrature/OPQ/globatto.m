% GLOBATTO Generalized Gauss-Lobatto quadrature formula
%
%    Given a measure by its (N+2r)x2 array ab of recurrence coefficients and the 
%    fixed nodes endl and endr, both of multiplicity r, this routine generates
%    the generalized Gauss-Lobatto quadrature rule involving r-1 consecutive 
%    derivative values at the fixed points endl and endr, and N interior 
%    points. The N interior nodes and weights are stored in the last N rows
%    of the (N+2r)x2 output array xw, and the fixed nodes endl and endr
%    (repeated r times) and the corresponding 2r weights associated with the
%    r derivative values in the first 2r rows of xw. It is assumed that N>2 and r>0.
%
function xw=globatto(N,ab,r,endl,endr)
if r<1, error('r too small'), end
if N<3, error('N too small'), end
odd=2*floor(r/2)~=r;
%
% internal nodes and weights
%
ab1=ab;
for s=2*r:-1:r+1
  ab0=ab1;
  ab1=chri1(N+s-1,ab0,endl);
end
for s=r:-1:1
  ab0=ab1;
  ab1=chri1(N+s-1,ab0,endr);
end
ab1(1,2)=abs(ab1(1,2)); xw0=gauss(N,ab1);
xw(2*r+1:N+2*r,:)=xw0;
xw(2*r+1:N+2*r,2)=xw(2*r+1:N+2*r,2)./(((xw(2*r+1:N+2*r,1)...
  -endl).*(endr-xw(2*r+1:N+2*r,1))).^r);
xw(1:r,1)=endl*ones(r,1);
xw(r+1:2*r,1)=endr*ones(r,1);
%
% external weights
%
r0=floor((r+1)/2); rs=zeros(r-1,1); tau=zeros(N+r0,1);
tau(1:N)=xw0(:,1);
for e=1:2
  if e==1
    a=endl; b=endr; tau(N+1:N+r0)=b;
  else
    a=endr; b=endl; tau(N+1:N+r0)=b;
  end
  A=zeros(r); rhs=zeros(r,1); x=zeros(r,1);
  for rho=1:r-1
    rs(rho)=sum((tau-a).^(-rho));
    if odd, rs(rho)=rs(rho)-.5*(b-a)^(-rho); end
  end
  for i=r:-1:1
    A(i,i)=factorial(i-1)*prod((tau-a).^2);
    if odd, A(i,i)=A(i,i)/(b-a); end
    if i<r
      for j=i+1:r
        A(i,j)=-2*sum(A(i+1:j,j).*rs(1:j-i))/(j-i);
      end
    end
  end
  ng=N+r;
  xwg=gauss(ng,ab);
  for i=1:r
    s=0;
    for ig=1:ng
      t=xwg(ig,1); w=xwg(ig,2);
      p=prod((tau-t).^2);
      if odd, p=p/(b-t); end 
      s=s+w*((t-a)^(i-1))*p;
    end
    rhs(i)=s;
  end
  x(r)=rhs(r)/A(r,r);
  for i=r-1:-1:1
    x(i)=(rhs(i)-sum(A(i,i+1:r).*(x(i+1:r))'))/A(i,i);
  end
  xw(1+(e-1)*r:e*r,2)=x;
end
