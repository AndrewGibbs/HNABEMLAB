% CAUCHY Cauchy integral of orthogonal polynomials.
%
%    Given a weight function w(t) through the first numax recurrence
%    coefficients ab of its orthogonal polynomials, [rho,r,nu]=CAUCHY
%    (N,ab,z,eps0,nu0,numax) generates to a relative accuracy eps0
%    the Cauchy integrals of the first N+1 orthogonal polynomials,
%    i.e., the N+1 integrals 
%
%        rho_k(z,w)=int pi_k(t)w(t)dt/(z-t),   k=0,1,...,N,
%
%    and their ratios
%
%              r_k=rho_{k+1}/rho_k,   k=0,1,...,N
%
%    where r_0=rho(0) and pi_k is the monic polynomial of degree k
%    orthogonal relative to the weight function w. The procedure 
%    used is a continued fraction algorithm, which involves backward 
%    recursion from nu down to 0. The starting index nu is increased
%    by 5 until convergence of the continued fraction algorithm is 
%    obtained. The initial value of nu to be used is nu0, and nu
%    is not allowed to exceed numax. The output variable nu is 
%    either numax, if there is no convergence, or the final value
%    of nu achieving convergence. The alpha- and beta-coefficients 
%    of the given weight function are to be provided in the first 
%    and second column of the (numax)x2 input array ab; the Cauchy 
%    integrals are stored in the (N+1)x1 output array rho and their
%    ratios in the (N+1)x1 output array r.
%
function [rho,r,nu]=cauchy(N,ab,z,eps0,nu0,numax)
if N>=numax, error('input variable numax too small'), end
N0=size(ab,1); if N0<numax, error('input array ab too short'), end
if nu0>numax, error(' input variables nu0 and numax incompatible'), end
if nu0<N+1, nu0=N+1; end
nu=nu0-5;
%nu=nu0-1;
rold=ones(N+1,1);
rho=zeros(N+1,1);
r=zeros(N+1,1);
while any(abs(r-rold)>eps0*abs(r))
  nu=nu+5;
%  nu=nu+1;
  if nu>numax, error('backward recurrence index exceeds limit'), end
  rold=r;
  rnu=0;
  for n=nu:-1:1
    rnu=ab(n,2)/(z-ab(n,1)-rnu);
    if n<=N+1, r(n)=rnu; end
  end
end
rho(1)=r(1);
if N==0, return, end
for n=2:N+1, rho(n)=r(n)*rho(n-1); end
