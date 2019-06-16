% KRONROD_RATIONAL Rational Gauss-Kronrod quadrature rule
%
%    The call xw=KRONROD_RATIONAL(N,abmod) generates from the first
%    2N+1  modified recurrence coefficients the (2N+1)-point rational
%    Gauss-Kronrod formula xw according to Theorem 3.41.
%
function xw=kronrod_rational(N,abmod)
global Z M
xw=kronrod(N,abmod);
if M==0, return, end
for n=1:2*N+1
  p(n)=prod((1+xw(n,1)*Z(1:M,1)).^Z(1:M,2));
end
xw(:,2)=xw(:,2).*p';
