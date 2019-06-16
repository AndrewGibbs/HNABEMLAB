% R_MULTIDOMAIN_STI A multidomain Stieltjes algorithm.
%
%    Given a finite number d of intervals (disjoint or not) and
%    on each interval a measure and corresponding Jacobi matrix,
%    the call ab=r_multidomain_sti(N,abmd) generates the first N
%    recurrence coefficients for the polynomials orthogonal on
%    the union of these intervals and measures. The d given
%    measures are identified via the Nx(2*d) input array abmd
%    holding the respective alpha- and beta-coefficients in two
%    consecutive columns of abmd. The desired recursion coefficients
%    are stored in the Nx2 output array ab, the first column
%    containing the alpha-, the second column the beta-coefficients.
%    The routine uses a matrix formulation of the Stieltjes
%    algorithm in combination with Gauss quadratures.
%
function ab=r_multidomain_sti(N,abmd)
d=size(abmd,2)/2;
J=zeros(N,d*N); I=eye(N);
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
Nd=1:N:d*N;
z=zeros(d*N,1); z0=z; z(Nd)=1;
num=0; den=0;
for j=1:d
  Nj=1+(j-1)*N:j*N;
  num=num+abmd(1,2*j)*(z(Nj)'*J(:,Nj)*z(Nj));
  den=den+abmd(1,2*j);
end
ab(1,1)=num/den; ab(1,2)=den;
for k=1:N-1
  zm1=z0; z0=z;
  for j=1:d
    Nj=1+(j-1)*N:j*N;
    z(Nj)=(J(:,Nj)-ab(k,1)*I)*z0(Nj)-ab(k,2)*zm1(Nj);
  end
  numa=0; numb=0;
  for j=1:d
    Nj=1+(j-1)*N:j*N;
    numb=numb+abmd(1,2*j)*(z(Nj)'*z(Nj));
    numa=numa+abmd(1,2*j)*(z(Nj)'*J(:,Nj)*z(Nj));
  end
  ab(k+1,1)=numa/numb; ab(k+1,2)=numb/den;
  den=numb;
end
