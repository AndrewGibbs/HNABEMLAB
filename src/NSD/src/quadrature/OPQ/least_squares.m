% LEAST_SQUARES Polynomial least squares approximation.
%
%    This routine generates the Nx(n+1) array phat of least squares
%    approximants evaluated at the abscissae of the discrete inner
%    product, and the (n+1)-vector c of Fourier coefficients. The
%    abscissae and weights of the inner product are input via the
%    first, resp. second, column of the Nx2 array xw. The (n+1)x2
%    input array ab is to contain the first n+1 alpha- and 
%    beta-coefficients in the recurrence relation for the MONIC
%    discrete orthogonal polynomials. The actual discrete orthogonal
%    polynomials used are normalized to have leading coefficients
%    given by the input vector d of dimension n+1, that is, the
%    coefficient of the leading power k is d(k+1), k=0,1,...,n.
%    The data values are given in the column vector f of dimension N.
%
function [phat,c]=least_squares(n,f,xw,ab,d)
N=size(xw,1);
if n>N-1, error('n too large'), end
%
% Generate the matrix of orthogonal polynomials
%
p=zeros(N,n+1); p2=zeros(1,n+1); pn=ab(1,2);
p(:,1)=d(1); p2(1)=d(1)^2*pn; pn=ab(2,2)*pn;
p(:,2)=d(2)*(xw(:,1)-ab(1,1)); p2(2)=d(2)^2*pn;
for k=2:n
  pn=ab(k+1,2)*pn;
  p(:,k+1)=d(k+1)*((xw(:,1)-ab(k,1)).*p(:,k)/d(k)-ab(k,2)...
  *p(:,k-1)/d(k-1));
  p2(k+1)=d(k+1)^2*pn;
end
%
% Compute the matrix of least squares approximants and the array
% of Fourier coefficients
%
phat=zeros(N,n+1); e=f;
c(1)=sum(xw(:,2).*e.*p(:,1))/p2(1); phat(:,1)=c(1)*p(:,1);
for i=2:n+1
  c(i)=sum(xw(:,2).*e.*p(:,i))/p2(i);
  phat(:,i)=phat(:,i-1)+c(i)*p(:,i);
  if i==n+1, return, end
  e=e-c(i)*phat(:,i);
end
