% G_N The function g_n.
%
%    This function is used to compute the absolute condition number
%    for the map G_n from the first 2n normalized modified moments
%    of a measure dlambda to the n-point Gauss quadrature formula.
%    The nodes and weights of the Gauss formula are input through
%    the nx2 array xw. The input variable t must be a
%    single variable; it cannot be an array.
%
function g=g_n(n,t,xw)
r=1:n; rho=zeros(n,1); sig=zeros(n,1);
for nu=1:n
  mu=find(r-nu);
  rho(nu)=prod(xw(nu,1)-xw(mu,1));
  sig(nu)=sum(1./(xw(nu,1)-xw(mu,1)));
end
g=sum((1./((t-xw(:,1)).*rho)).^4.*((1-2*(t-xw(:,1))...
  .*sig).^2+((t-xw(:,1))./xw(:,2)).^2))/sum(1./((t-...
  xw(:,1)).*rho))^4;
