% SPJAC Stieltjes-Perron inversion for Jacobi polynomials
%
%    This routine tests the Stieltjes-Perron inversion formula,
%    as implemented by the routine SPHT, on Jacobi measures with
%    parameters alpha and beta in (-1,1], recovering the Jacobi
%    weight function for selected values of x in (-1,1) from the
%    recurrence coefficients.
%
f1='                        x=%4.2f\n';
f2='     %4.0f  %4.0f  %3.1f  %4.1f %10.3e %10.3e\n';
n=9; numax=4000; iopt=1; eps0=1e-6; err=zeros(n,n+1);
for x=[-.99 -.9 -.8 -.4 0 .4 .8 .9 .99]
  fprintf(f1,x)
  fprintf('\n')
  disp('      nu0   nu   alpha beta   err(1,9)  err(10,1)')
  for alpha=-.8:.2:1
    for beta=alpha:.2:1
      ab=r_jacobi(numax,alpha,beta);
      nu0=3880;
      [E,w,nu0,nu]=SPHT(n,x,iopt,ab,nu0,numax,eps0);
      exact=(1-x)^alpha*(1+x)^beta;
      err(n,2)=abs(E(n,2)-exact);
      err(1,n+1)=abs(E(1,n+1)-exact);
      fprintf(f2,nu0,nu,alpha,beta,err(n,2),err(1,n+1))
    end
  end
end
