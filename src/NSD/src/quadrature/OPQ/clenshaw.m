% CLENSHAW Clenshaw algorithm.
%
%    Given functions y_k(x), k=0,1,...,N, N>=0, satisfying a three-term
%    recurrence relation
%
%  y_{k+1}(x)=(x-alpha_k) y_k(x) - beta_k y_{k-1}, k=0,1,...,N-1,
%
%    with initial values y_0, y_{-1}, the call s=CLENSHAW(N,x,y0,ym1,ab,c) uses
%    Clenshaw's algorithm to evaluate the sum
%
%                   s=sum_{j=0}^N c_j y_j(x).
%
%    The initial values y_0, y_{-1} are input through the parameters y0, ym1, 
%    and the alpha- and beta-coefficients (if N>0) through the Nx2 array ab, 
%    the first column of which contains the alpha's, and the second column 
%    the beta's. The coefficients c_{k-1} are input via the 1x(N+1) array c.
%    The argument x is allowed to be vector-valued.
%
function s=clenshaw(N,x,y0,ym1,ab,c)
if N<0, error('N out of range'), end
if size(x,2)>1, x=x'; end
lx=length(x);
if N==0, s=c(1)*y0*ones(lx,1); return, end
if size(ab,1)<N, error('array ab too short'), end
if size(c,2)<N+1, error('array c too short'), end
u=zeros(lx,N+1);   
u(:,N+1)=c(N+1); u(:,N)=(x-ab(N,1)).*u(:,N+1)+c(N);
for k=N-1:-1:1
  u(:,k)=(x-ab(k,1)).*u(:,k+1)-ab(k+1,2).*u(:,k+2)+c(k);
end
s=u(:,1)*y0-ab(1,2)*u(:,2)*ym1;
