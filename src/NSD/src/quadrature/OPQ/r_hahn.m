% R_HAHN Recurrence coefficients for monic Hahn polynomials.
%
%    ab=R_HAHN(N,a,b) generates the N+1 recurrence coefficients 
%    for the monic Hahn polynomials with parameters a and b. These
%    are orthogonal on the discrete set of N+1 points 0,1,2,...,N
%    with weights (a+k choose k)(b+N-k choose N-k), k=0,1,2,...,N.
%    The N+1 alpha-coefficients are stored in the first column,
%    the N+1 beta-coefficients in the second column, of the
%    (N+1)x2 array ab. The call ab=R_HAHN(N,a) is the same as ab=
%    R_HAHN(N,a,a), and ab=HAHN(N) the same as ab=R_HAHN(N,0,0),
%    which produces the recurrence coefficients of the discrete 
%    Chebyshev polynomials for the points 0,1,2,...,N.
%
function ab=r_hahn(N,a,b)
if nargin<2, a=0; end; if nargin<3, b=a; end
if((N<=0)|(a<-1)|(b<-1)), error('parameter(s) out of range'), end
n=1:N; ab(1,2)=prod(1+(a+b+1)./n);
if a+b==0
  n=(1:N+1)';
  ab(n,1)=((2*n+a+b-1)*N+(b-a)*n+a)./(2*(2*n-1));
  n=(1:N)';
  ab(n+1,2)=.25*((N+1)^2)*(1+a./n).*(1+b./n).*(1-(n./(N+1)).^2)...
    ./(4-(1./n).^2);
elseif a+b+1==0
  n=(1:N+1)';
  ab(n,1)=((2*(n-1).^2+b)*N+(2*b+1)*(n-1).^2)./(4*(n-1).^2-1);
  n=(1:N)';
  ab(n+1,2)=.25*((N+1)^2)*(1+a./n).*(1+b./n).*(1-n./(N+1)).*...
    (1+(n-1)./(N+1))./(4-(1./n).^2);
else
  n=(1:N+1)';
  ab(n,1)=((n+a+b).*(n+a).*(N-n+1)./(2*n+a+b)+(n-1).*(n+b-1).*...
    (N+n+a+b)./(2*n+a+b-2))./(2*n+a+b-1);
  n=(1:N)';
  ab(n+1,2)=((N+1)^2)*(1+a./n).*(1+b./n).*(1+(a+b)./n).*...
    (1-n./(N+1)).*(1+(n+a+b)./(N+1))./(((2+(a+b)./n).^2).*...
    ((2+(a+b)./n).^2-(1./n).^2));
end
