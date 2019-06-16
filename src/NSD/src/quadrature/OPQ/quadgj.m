% QUADGJ A quadrature routine used in R_GJACOBI.
%
function xw=quadgj(N,i)
global a b mc tc
ab=r_jacobi(N,tc(i+1,2),tc(i,2));
xw=gauss(N,ab);
xw(:,1)=.5*(tc(i+1,1)-tc(i,1))*xw(:,1)+.5*(tc(i+1,1)...
  +tc(i,1));
xw(:,2)=((.5*(tc(i+1,1)-tc(i,1)))^(tc(i+1,2)+...
  tc(i,2)+1))*xw(:,2);
j=1:mc+1;
nu=find((j-i).*(j-i-1));
for n=1:N
  xw(n,2)=xw(n,2)*prod(abs((xw(n,1)-tc(nu,1))).^tc(nu,2));
end
%
% If the generalized Jacobi weight function contains a
% function phi, then the following statement must be
% added
%
%  xw(:,2)=phi(xw(:,1)).*xw(:,2);
%
% The routine for the function phi must be written in
% such a way that it accepts vector-valued arguments.
%
