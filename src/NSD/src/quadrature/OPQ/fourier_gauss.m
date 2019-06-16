% FOURIER_GAUSS Fourier coefficients based on Gauss 
%    quadrature points and weights.
%
function c=fourier_gauss(N,f)
k=1:N; tk=cos((2.*k-1)*pi/(2*N)); 
Tm1=zeros(1,N); T0=ones(1,N); Tp1=tk;
c(1)=2*sum(f'.*T0)/N; c(2)=2*sum(f'.*Tp1)/N;
for n=2:N-1
  Tm1=T0; T0=Tp1;
  Tp1=2*tk.*T0-Tm1;
  c(n+1)=2*sum(f'.*Tp1)/N;
end
