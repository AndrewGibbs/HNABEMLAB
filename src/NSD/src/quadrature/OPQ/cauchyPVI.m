% CAUCHYPVI Cauchy principal value integral.
%
%    The call cpvi=cauchyPVI(N,x,f,ddf,iopt,ab,rho0) computes an 
%    N-point (resp. (N+1)-point, if iopt~=1) approximation to the 
%    Cauchy principal value integral int f(t)dlambda(t)/(x-t). 
%    The integration measure is input via the (N+1)x2 array ab
%    of the recurrence coefficients for the monic orthogonal
%    polynomials associated with dlambda, the alpha-coefficients 
%    being stored in the first column, the beta-coefficients in 
%    the second column of ab. The Hilbert transform of the 
%    measure, int dlambda(t)/(x-t), has to be input through the
%    parameter rho0. The routine assumes the function f to be
%    provided by the M-file f.m, and it also relies on the M-file 
%    ddf.m for providing a carefully programmed routine for the 
%    divided difference (f(x)-f(y))/(x-y) of f. Both routines 
%    must be able to accept vector-valued arguments.
%
function cpvi=cauchyPVI(N,x,f,ddf,iopt,ab,rho0)
if iopt==1
  N0=N;
else
  N0=N+1;
end
a=zeros(1,N0);
xw=gauss(N,ab);
PI=zeros(N+1,N); pi2=zeros(1,N+1); PID=zeros(1,N);
PI(1,:)=1; PI(2,:)=xw(:,1)'-ab(1,1);
pi2(1)=ab(1,2); pi2(2)=ab(2,2)*pi2(1);
if N>1
  for k=2:N
    PI(k+1,:)=(xw(:,1)'-ab(k,1)).*PI(k,:)-ab(k,2).*PI(k-1,:);
    pi2(k+1)=ab(k+1,2)*pi2(k);
  end
end
for k=1:N
  a(k)=sum(xw(:,2)'.*PI(k,:).*feval(f,xw(:,1)'))/pi2(k);
end
if iopt~=1
  for nu=1:N
    d=xw(nu,1)-xw(:,1); PID(nu)=prod(d(find(d)));
  end
  a(N0)=sum(feval(ddf,x*ones(1,N),xw(:,1)')./PID);
end
cpvi=clenshaw(N0-1,x,rho0,1,ab,a);
