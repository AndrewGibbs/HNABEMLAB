% FEJER  Fejer quadrature rule.
%        The call uv=fejer(n) generates the n-point Fejer quadrature rule.
function uv=fejer(N)
n=N:-1:1; m=1:floor(N/2);
th=(2*n-1)*pi./(2*N); 
uv(:,1)=cos(th');
for k=N:-1:1
  s=sum(cos(2*m*th(k))./(4*(m.^2)-1));
  uv(k,2)=2*(1-2*s)/N;
end
