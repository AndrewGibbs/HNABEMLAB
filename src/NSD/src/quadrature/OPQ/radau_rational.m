% RADAU_RATIONAL Rational Gauss-Radau quadrature rule.
%
%    The call xw=RADAU_RATIONAL(N,abmod,end0) generates from the first N+1 
%    modified recurrence coefficients the (N+1)-point rational Gauss-Radau 
%    quadrature xw according to the Theorem 3.39.
%
function xw=radau_rational(N,abmod,end0)
global Z M
xw=radau(N,abmod,end0);
if M==0, return, end
for n=1:N+1
  p(n)=prod((1+xw(n,1)*Z(1:M,1)).^Z(1:M,2));
end
xw(:,2)=xw(:,2).*p';

