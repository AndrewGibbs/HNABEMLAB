% KRONROD Gauss-Kronrod quadrature formula.
%
%    xw=KRONROD(n,ab) generates the (2n+1)-point Gauss-Kronrod
%    quadrature rule for the weight function w encoded by the 
%    recurrence matrix ab of order [ceil(3*n/2)+1]x2 containing
%    in its first and second column respectively the alpha- and 
%    beta-coefficients in the three-term recurrence relation
%    for w. The 2n+1 nodes, in increasing order, are output
%    into the first column, the corresponding weights into the
%    second column, of the (2n+1)x2 array xw.
%
%    Supplied by Dirk Laurie, 6-22-1998
%
function xw=kronrod(N,ab)
ab0=r_kronrod(N,ab);
if(sum((ab0(:,2)>0))~=2*N+1) error('Gauss-Kronrod does not exist'), end
J=zeros(2*N+1);
for k=1:2*N
  J(k,k)=ab0(k,1);
  J(k,k+1)=sqrt(ab0(k+1,2));
  J(k+1,k)=J(k,k+1);
end
J(2*N+1,2*N+1)=ab0(2*N+1,1);
[V,D]=eig(J);
d=diag(D);
e=ab0(1,2).*(V(1,:).^2);
[x,i]=sort(d);
w=e(i)';
xw=[x w];
