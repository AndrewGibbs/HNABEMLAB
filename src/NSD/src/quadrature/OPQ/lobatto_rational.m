% LOBATTO_RATIONAL Rational Gauss-Lobatto quadrature rule.
%
%    The call xw=LOBATTO_RATIONAL(N,abmod) generates from the first N+2 
%    modified recurrence coefficients the (N+2)-point rational Gauss-type 
%    quadrature xw according to the Theorem 3.40.
%
function xw=gauss_rational(N,abmod)
global Z M
xw=gauss(N,abmod);
if M==0, return, end
for n=1:N
  p(n)=prod((1+xw(n,1)*Z(1:M,1)).^Z(1:M,2));
end
xw(:,2)=xw(:,2).*p';

