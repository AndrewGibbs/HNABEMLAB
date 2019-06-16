% HTJAC Hilbert transform of the Jacobi measure
%
%    This computes the Hilbert transform of the Chebyshev measures
%    of the first and second kind for selected values of x from
%    the respective Cauchy integrals using a procedure implemented
%    in the routine SPHT.
%
f1='                        x=%4.2f\n';
f2='     %4.0f  %4.0f  %3.1f  %4.1f %10.3e %10.3e\n';
n=9; iopt=2; numax=4000; eps0=1e-6; err=zeros(n,n+1);
for x=[0 .2 .4 .6 .8 .9 .99]
  fprintf(f1,x)
  fprintf('\n')
  disp('      nu0   nu   alpha beta   err(1,9)  err(10,1)')
  for alpha=[-.5 .5]
  beta=alpha;
    ab=r_jacobi(numax,alpha,beta);
    nu0=3880;
    [E,w,nu0,nu]=SPHT(n,x,iopt,ab,nu0,numax,eps0);
    if alpha==-.5
      exact=0;
    else
      exact=-pi*x;
    end
    for k=2:2:n+1
      for m=1:n+2-k
        err(m,k)=abs(E(m,k)-exact);
      end
    end
    fprintf(f2,nu0,nu,alpha,beta,err(n,2),err(1,10))
  end
end
