% GRADAU Generalized Gauss-Radau quadrature formula
%
%    Given a measure by its (N+r)x2 array of recurrence 
%    coefficients and the fixed node endl of multiplicity r,
%    this routine generates the genralized Gauss-Radau
%    quadrature rule involving r-1 consecutive derivative 
%    values at the fixed point endl and N interior points.
%    The N interior nodes and weights are stored in the
%    last N rows of the (N+r)x2 output array xw, and the
%    fixed node endl (repeated r times) and the r weights 
%    associated with the r derivative values in the first
%    r rows of xw.
%
function xw=gradau(N,ab,r,endl)
%
% internal nodes and weights
%
ab1=ab;
for s=r:-1:1
  ab0=ab1;
  ab1=chri1(N+s-1,ab0,endl);
end
xw0=gauss(N,ab1);
xw(r+1:N+r,:)=xw0;
xw(r+1:N+r,2)=xw(r+1:N+r,2)./((xw(r+1:N+r,1)-endl).^r);
xw(1:r,1)=endl*ones(r,1);
%
% external weights
%
rs=zeros(r-1,1); 
A=zeros(r); b=zeros(r,1); x=zeros(r,1);
for rho=1:r-1
  rs(rho)=sum((xw0(:,1)-endl).^(-rho));
end
for i=r:-1:1
  A(i,i)=factorial(i-1)*prod((endl-xw0(:,1)).^2);
  if i<r
    for j=i+1:r
      A(i,j)=-2*sum(A(i+1:j,j).*rs(1:j-i))/(j-i);
    end
  end
end
ng=N+floor((r+1)/2);
xwg=gauss(ng,ab);
for i=1:r
  s=0;
  for ig=1:ng
    t=xwg(ig,1); w=xwg(ig,2);
    s=s+w*((t-endl)^(i-1)*prod((t-xw0(:,1)).^2));
  end
  b(i)=s;
end
x(r)=b(r)/A(r,r);
for i=r-1:-1:1
  x(i)=(b(i)-sum(A(i,i+1:r).*(x(i+1:r))'))/A(i,i);
end
xw(1:r,2)=x;
