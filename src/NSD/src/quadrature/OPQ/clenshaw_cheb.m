% CLENSHAW_CHEB Clenshaw algorithm for Chebyshev series.
%
%    The call s=CLENSHAW_CHEB(N,x,c) uses Clenshaw's algorithm to evaluate 
%    the finite Chebyshev series
%
%           s=.5 c_0 + sum_{j=1}^N c_j T_j(x),  N>=0,
%
%    the coefficients c_{j-1} being input via the 1x(N+1) array c.
%    The argument x is allowed to be vector-valued.
%
function s=clenshaw_cheb(N,x,c)
if N<0, error('N out of range'), end
if size(x,2)>1, x=x'; end
lx=length(x);
if N==0, s=.5*c(1)*ones(lx,1); return, end
if size(c,2)<N+1, error('array c too short'), end
v=zeros(lx,N+2);
v(:,N+1)=c(N+1);
for k=N:-1:1
  v(:,k)=2*x.*v(:,k+1)-v(:,k+2)+c(k);
end
s=.5*(v(:,1)-v(:,3));
