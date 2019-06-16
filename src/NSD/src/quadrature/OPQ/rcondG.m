% RCONDG Relative condition number of G_n.
%
%    This computes the relative condition number of the map G_n
%    from the first 2n moments of a measure dlambda to the n-point
%    Gauss quadrature formula. The measure is specified via the
%    first n recursion coefficients, stored in the nx2 array ab,
%    of the polynomials orthogonal with respect to the measure. If
%    iopt=1, the moments, input through the (2n)-vector mom, are
%    ordinary moments, otherwise modified moments relative to
%    orthogonal polynomials specified via their first 2n recursion
%    coefficients, stored in the (2n)x2 input array ab0. If iopt=1,
%    use the zero array for ab0.
%
function rc=rcondG(n,iopt,ab,ab0,mom)
xw=gauss(n,ab); C=zeros(2*n);
if iopt==1
  U=Tinv(n,xw(:,1));
  for kap=1:n
    for lam=1:2*n
      C(kap,lam)=moll(mom(lam))*U(kap,lam)/moll(xw(kap,2));
      C(n+kap,lam)=moll(mom(lam))*U(n+kap,lam)/...
        (xw(kap,2)*moll(xw(kap,1)));
    end
  end
  rc=norm(C,inf);
else
  xw0=gauss(2*n,ab0);
  rho=zeros(n,1); sig=zeros(n,1);
  for nu=1:n
    d=xw(nu,1)-xw(:,1);
    rho(nu)=prod(d(find(d))); sig(nu)=1/sum(d(find(d)));
  end
  for k=1:2*n
    pnorm2(k)=prod(ab0(1:k,2));
  end
  momnorm=mom(1:2*n)./sqrt(pnorm2);
  for nu=1:n
    for mu=1:2*n
      t=xw0(:,1);
      p0=zeros(2*n,1); p=ones(2*n,1);
      if mu>1
        for k=1:mu-1
          pm1=p0; p0=p;
          p=(t-ab0(k,1)).*p0-ab0(k,2).*pm1;
        end
      end
      den=zeros(2*n,1);
      for kap=1:n
        den=den+((t-xw(kap,1))*rho(kap)).^(-1);
      end
      den=den.^2;
      hherm=(1-2*(t-xw(nu,1))*sig(nu)).*((t-xw(nu,1))...
        *rho(nu)).^(-2)./den;
      kherm=((t-xw(nu,1))*rho(nu)).^(-2).*(t-xw(nu,1))./den;
      a=sum(xw0(:,2).*p.*hherm)/pnorm2(mu);
      b=sum(xw0(:,2).*p.*kherm)/pnorm2(mu);
      C(nu,mu)=sqrt(pnorm2(mu))*moll(momnorm(mu))*abs(a)/...
        moll(xw(nu,2));
      C(n+nu,mu)=sqrt(pnorm2(mu))*moll(momnorm(mu))*abs(b)/...
        (xw(nu,2)*moll(xw(nu,1)));
    end
  end
  rc=norm(C);
end
