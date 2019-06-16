% CLENSHAW_SOB Clenshaw-like algorithm for Sobolev orthogonal polynomials
%
%    Given the NxN matrix B=B_N of the Sobolev recurrence
%    coefficients, s=CLENSHAW_SOB(n,N,s,x,B,c) generates the sum
%    r=c_0 pi_0(x)+c_1 pi_1(x)+...+c_n pi_n(x) in (monic) Sobolev 
%    orthogonal polynomials and its first s derivatives, where
%    s<=2. The coefficients c_j are assumed input via the 1x(n+1) 
%    array c. The argument x is allowed to be vector-valued.
%
function r=clenshaw_sob(n,N,s,x,B,c)
if N<0|n<0|n>N, error('N or n out of range'), end
if s>2, error('s too large'), end
if size(x,2)>1, x=x'; end
lx=length(x);
if max(size(B)<[N N])>0, error('array B too small'), end
if size(c,2)<n+1, error('array c too short'), end
r=zeros(lx,s+1);
u=zeros(lx,n+1); v=zeros(lx,n+1); w=zeros(lx,n+1);
u(:,n+1)=c(n+1); v(:,n+1)=0; w(:,n+1)=0;
for k=n:-1:1
  u(:,k)=(x-B(1,k)).*u(:,k+1)+c(k);
  if k<n
    for l=1:n-k
      u(:,k)=u(:,k)-B(l+1,k+l).*u(:,k+l+1);
    end
  end
end
r(:,1)=u(:,1);
if s>0
  for k=n:-1:1
    v(:,k)=(x-B(1,k)).*v(:,k+1)+u(:,k+1);
    if k<n 
      for l=1:n-k
        v(:,k)=v(:,k)-B(l+1,k+l).*v(:,k+l+1);
      end
    end
  end
  r(:,2)=v(:,1);
  if s>1
    for k=n:-1:1
      w(:,k)=(x-B(1,k)).*w(:,k+1);
      if k>1, w(:,k)=w(:,k)+u(:,k+1); end
      if k<n
        for l=1:n-k
          w(:,k)=w(:,k)-B(l+1,k+l).*w(:,k+l+1);
        end
      end
    end
    for k=n:-1:1
      v(:,k)=(x-B(1,k)).*v(:,k+1)+w(:,k+1);
      if k<n
        for l=1:n-k
          v(:,k)=v(:,k)-B(l+1,k+l).*v(:,k+l+1);
        end
      end
    end
    r(:,3)=2*v(:,1);
  end
end
