% FOURIER_LOBATTO Fourier coefficients based on Gauss-Lobatto 
%    quadrature points and weights.
%
function c=fourier_lobatto(N,f)
k=1:N; tk=cos((N-k)*pi/(N-1)); wk=ones(1,N);
wk(1)=.5; wk(N)=wk(1);
Tm1=zeros(1,N); T0=ones(1,N); Tp1=tk;
c(1)=2*sum(wk.*f'.*T0)/(N-1); 
c(2)=2*sum(wk.*f'.*Tp1)/(N-1);
for n=2:N-1
  Tm1=T0; T0=Tp1;
  Tp1=2*tk.*T0-Tm1;
  c(n+1)=2*sum(wk.*f'.*Tp1)/(N-1);
end
