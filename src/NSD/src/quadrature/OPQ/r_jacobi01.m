% R_JACOBI01 Recurrence coefficients for monic Jacobi polynomials 
% on [0,1].
%
%    ab=R_JACOBI01(n,a,b) generates the first n recurrence
%    coefficients for monic Jacobi polynomials on [0,1] with 
%    parameters a and b. These are orthogonal on [0,1] relative 
%    to the weight function w(t)=(1-t)^a t^b. The n alpha-
%    coefficients are stored in the first column, the n beta-
%    coefficients in the second column, of the nx2 array ab. The
%    call ab=R_JACOBI01(n,a) is the same as ab=R_JACOBI01(n,a,a)
%    and ab=R_JACOBI01(n) the same as ab=R_JACOBI01(n,0,0).
%
function ab=r_jacobi01(N,a,b)
if nargin<2, a=0; end;  if nargin<3, b=a; end
if((N<=0)|(a<=-1)|(b<=-1)) error('parameter(s) out of range'), end
cd=r_jacobi(N,a,b);
n=1:N;  ab(n,1)=(1+cd(n,1))./2;  
ab(1,2)=cd(1,2)/2^(a+b+1);
n=2:N;  ab(n,2)=cd(n,2)./4;
