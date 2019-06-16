% QUADJP A quadrature routine used in R_JACPLUS.
%
function xw=quadjp(N,i)
global a b
ab=r_jacobi(N,a,b); ab(1,2)=1;
xw=gauss(N,ab);
