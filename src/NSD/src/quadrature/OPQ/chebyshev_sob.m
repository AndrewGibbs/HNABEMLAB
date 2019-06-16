% CHEBYSHEV_SOB Modified Chebyshev algorithm for Sobolev orthogonal polynomials
%
%    For given n>0, let the 2x(2*n) array mom contain the first 2*n
%    modified moments of the two measures dlambda_0 and dlambda_1,
%    and the (2*n-1)x2 array abm the first 2*n-1 recurrence
%    coefficients of the orthogonal polynomials defining the
%    modified moments. The call [B,normsq]=chebyshev_sob(n,mom,abm)
%    generates the nxn matrix B of the recurrence coefficients
%    of the monic Sobolev orthogonal polynomials (with s=1), with
%    beta_j^k, 0<=j<=k, k=0,1,...,n-1, occupying the position
%    B(j+1,k+1) of the matrix B and all remaining elements of B
%    being zero. The routine also returns the n-vector normsq of
%    the squared norms of the Sobolev orthogonal polynomials.
%    The call [B,normsq]=chebyshev_sob(n,mom] does the same, but 
%    using ordinary moments.
%
%    REFERENCE: Section 2 of W. Gautschi and M. Zhang, ``Computing 
%    orthogonal polynomials in Sobolev spaces'', Numer. Math. 71
%    (1995), 159-183.
%
function [B,normsq]=chebyshev_sob(n,mom,abm)
if n<=0, error('n out of range'), end
if n>size(mom,2)/2, n=size(mom,2)/2; end
if nargin<3, abm=zeros(2*n-1,2); end
if n>(size(abm,1)+1)/2, n=(size(abm,1)+1)/2; end
B=zeros(n); sig=zeros(n,2*n);
% 
% Initialization
%
sig(1,:)=mom(1,:); B(1,1)=abm(1,1)+sig(1,2)/sig(1,1);
normsq(1)=sig(1,1);
if n==1, return, end
mu(1,1,1,:)=mom(2,:);
mu(1,2,1,1)=0; mu(2,1,1,1)=0; mu(2,1,1,2)=0; mu(1,2,1,2)=mom(2,1);
tau=zeros(3,2*n-1); tau(3,1)=1;
for k=2:2*n-2
  for l=1:k-1
    tau(1,l)=tau(2,l); tau(2,l)=tau(3,l);
  end
  s=0;
  tau(3,1)=tau(2,2)*abm(2,2)+tau(2,1)*(abm(1,1)-abm(k,1))...
    -tau(1,1)*abm(k,2);
  s=s+tau(3,1)*mom(2,1);
  for j=2:k
    tau(3,j)=tau(2,j-1)+tau(2,j+1)*abm(j+1,2)+tau(2,j)*...
      (abm(j,1)-abm(k,1))-tau(1,j)*abm(k,2);
    if j==k, tau(3,j)=tau(3,j)+1; end
    s=s+tau(3,j)*mom(2,j);
  end
  mu(1,2,1,k+1)=s; mu(2,1,1,k+1)=0;
end
%
% Continuation
%
for k=2:n
  for l=k:2*n-k+1
    x=B(1:k-1,k-1); y=sig(k-1:-1:1,l);
    s=sum(x.*y);
    sig(k,l)=sig(k-1,l+1)+abm(l,1)*sig(k-1,l)...
      +abm(l,2)*sig(k-1,l-1)+mu(1,2,k-1,l)...
      -mu(2,1,k-1,l)-s;
    if(k==l), normsq(k)=sig(k,l); end
  end
  for j=1:k
    B(k-j+1,k)=(sig(j,k+1)+abm(k,1)*sig(j,k))/sig(j,j);
    if j<k
      B(k-j+1,k)=B(k-j+1,k)+abm(k,2)*sig(j,k-1)/sig(j,j);
      s=0;
      for l=j:k-1
        s=s+B(l-j+1,l)*sig(l,k)/sig(l,l);
      end
      B(k-j+1,k)=B(k-j+1,k)-s;
    end
    if j>1
      B(k-j+1,k)=B(k-j+1,k)-sig(j-1,k)/sig(j-1,j-1);
    end
  end
  if k<n
    for l=k:2*n-k
      for i2=1:2
        for i1=1:2
          if i1*i2<4
            s=0;
            for j=1:k-1
              s=s+B(j,k-1)*mu(i1,i2,k-j,l);
            end
            mu(i1,i2,k,l)=mu(i1,i2,k-1,l+1)+abm(l,1)*...
              mu(i1,i2,k-1,l)+abm(l,2)*mu(i1,i2,k-1,l-1)-s;
            if i1~=i2
              if i2==1
                sgn=1;
              else
                sgn=-1;
              end
              mu(i1,i2,k,l)=mu(i1,i2,k,l)+sgn*mu(1,1,k-1,l);
            end
          end
        end
      end
    end
  end
end

