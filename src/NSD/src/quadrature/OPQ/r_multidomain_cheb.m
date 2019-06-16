% R_MULTIDOMAIN_CHEB A multidomain modified Chebyshev algorithm.
%
%    Given a finite number d of intervals (disjoint or not) and
%    on each interval a measure and corresponding Jacobi matrix,
%    the call ab=r_multidomain_cheb(N,abmd,abmm) generates the 
%    first N recurrence coefficients for the polynomials orthogonal 
%    on the union of these intervals and measures. The d given
%    measures are identified via the Nx(2*d) input array abmd
%    holding the respective alpha- and beta-coefficients in two
%    consecutive columns of abmd. The desired recursion coefficients
%    are stored in the Nx2 output array ab, the first column
%    containing the alpha-, the second column the beta-coefficients.
%    The routine uses a matrix formulation of the modified Chebyshev
%    algorithm in combination with Gauss quadratures. The 2*N-1
%    recurrence coefficients of the polynomials defining the
%    modified moments are input through the (2*N-1)x2 array abmm.
%
function ab=r_multidomain_cheb(N,abmd,abmm)
d=size(abmd,2)/2;
J=zeros(N,d*N); I=eye(N); mom=zeros(1,2*N);
me=2:2:2*d;
mo=1:2:2*d-1;
for n=1:N
  m=n:N:n+(d-1)*N;
  J(n,m)=abmd(n,mo);
  if n<N
    J(n,m+1)=sqrt(abmd(n+1,me));
  end
  if n>1
    J(n,m-1)=J(n-1,m);
  end
end
mom(1)=sum(abmd(1,2:2:2*d));
Nd=1:N:d*N;
z=zeros(d*N,1); z0=z; z(Nd)=1; e1=z;
for k=1:2*N-1
  zm1=z0; z0=z;
  for j=1:d
    Nj=1+(j-1)*N:j*N;
    z(Nj)=(J(:,Nj)-abmm(k,1)*I)*z0(Nj)-abmm(k,2)*zm1(Nj);
  end
  s=0;
  for j=1:d
    Nj=1+(j-1)*N:j*N;
    s=s+abmd(1,2*j)*(z(Nj)'*e1(Nj));
  end
  mom(k+1)=s;
end
ab=chebyshev(N,mom,abmm);
