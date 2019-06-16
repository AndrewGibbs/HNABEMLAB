% TURAN Gauss-Turan quadrature formula.
%
%    The call xw=TURAN(n,s,eps0,ab0,hom) computes the n-point
%    Gauss-Turan quadrature formula for a measure dlambda, where
%    the quadrature sum has terms involving successive derivatives
%    of orders up to, and including, 2s. The input parameter
%    eps0 is an error tolerance. The measure dlambda is specified
%    by the ((s+1)*n)x2 array ab0 of the recurrence coefficients
%    for the orthogonal polynomials belonging to the measure.
%    If hom=1, a discrete homotopy in the variable s is used to
%    facilitate convergence of Newton's method for computing the
%    quadrature nodes. Otherwise, no homotopy is used. The nodes
%    of the quadrature rule are stored in the first column of the
%    nx(2*s+2) array xw, the successive weights in the remaining
%    2*s+1 columns.
%
function xw=turan(n,s,eps0,ab0,hom)
if hom==1, s0=1; else s0=s; end
xw=zeros(n,2*s+2);
sn=(s+1)*n; sn2=(s+1)*n^2;
if size(ab0,1)<sn, error('array ab0 too short'), end
xw0=gauss(sn,ab0);
f=zeros(2*n-1,1); fjac=zeros(2*n-1);
A=zeros(n+1,sn2); B=zeros(n+1,sn2); D=zeros(n,sn2);
P=zeros(n,sn2); Q=zeros(n,sn2);
PI=zeros(n+1,sn); PI(1,:)=ones(1,sn);
ab=ab0(1:n,:);
for sp=s0:s
  rho=[ab(:,1);ab(2:n,2)]; rhop=zeros(2*n-1,1);
%
% Newton iteration
%
  it=0;
  while any(abs(rho-rhop)>eps0*(abs(rho)+1))
    rhop=rho;
    it=it+1;
    if it>50, error('Newton iteration does not converge'), end
    abp=[rhop(1:n) [0;rhop(n+1:2*n-1)]];
    PI(2,:)=xw0(:,1)'-abp(1,1)*PI(1,:);
    for k=2:n
      PI(k+1,:)=(xw0(:,1)'-abp(k,1)).*PI(k,:)-abp(k,2)*PI(k-1,:);
    end
%
% Computation of f
%
    for nu=1:n
      f(2*nu-1)=((abp(nu,1)-xw0(:,1)').*PI(nu,:).^2.*PI(n+1,:).^(...
        2*sp))*xw0(:,2);
      if nu<n
        f(2*nu)=((abp(nu+1,2)*PI(nu,:).^2-PI(nu+1,:).^2).*PI(...
          n+1,:).^(2*sp))*xw0(:,2);
      end
    end  
%
% Computation of auxiliary matrices
%
    for mu=1:n
      m=mu:n:sn2;
      A(mu+1,m)=-PI(mu,:);
      for nu=mu+1:n
        A(nu+1,m)=(xw0(:,1)'-abp(nu,1)).*A(nu,m)-abp(nu,2).*A(...
          nu-1,m);
      end
      if mu>1
        B(mu+1,m)=-PI(mu-1,:);
        for nu=mu+1:n
          B(nu+1,m)=(xw0(:,1)'-abp(nu,1)).*B(nu,m)-abp(nu,2)*B(...
            nu-1,m);
        end
      end
    end
    for nu=1:n
      for mu=1:n
        m=mu:n:sn2;
        if nu==mu, D(nu,m)=1; end
        P(nu,m)=PI(nu,:).*(A(nu,m).*PI(n+1,:)+sp*A(n+1,m).*PI(...
          nu,:));
        Q(nu,m)=PI(nu,:).*(B(nu,m).*PI(n+1,:)+sp*B(n+1,m).*PI(...
          nu,:));
      end
    end
%
% Computation of the Jacobian of f
%
    for nu=1:n
      for mu=1:n
        m=mu:n:sn2;
        fjac(2*nu-1,mu)=2*(PI(n+1,:).^(2*sp-1).*((abp(nu,1)...
          -xw0(:,1)').*P(nu,m)+.5*D(nu,m).*PI(nu,:).^2.*PI(...
          n+1,:))*xw0(:,2));
        if mu>1
          fjac(2*nu-1,mu+n-1)=2*(PI(n+1,:).^(2*sp-1).*(abp(nu,1)...
            -xw0(:,1)').*Q(nu,m))*xw0(:,2);
        end
        if nu<n
          fjac(2*nu,mu)=2*(PI(n+1,:).^(2*sp-1).*(abp(nu+1,2)...
            *P(nu,m)-P(nu+1,m)))*xw0(:,2);
          if mu>1
            fjac(2*nu,mu+n-1)=2*(PI(n+1,:).^(2*sp-1).*((abp(nu+1,2)...
              *Q(nu,m)-Q(nu+1,m))+.5*D(nu+1,mu).*PI(nu,:).^2.*PI(...
                n+1,:))*xw0(:,2));
          end
        end
      end
    end
%
% Newton correction
%
    d=fjac\f; rho=rhop-d;
  end
  ab=[rho(1:n) [0;rho(n+1:2*n-1)]];
end
%it 
xwg=gauss(n,ab);
xw(:,1)=xwg(:,1);
%
% Computation of the weights
%
Ahat=zeros(2*s+1);
dtau=zeros(n,1); u=zeros(1,2*s);
for nu=1:n
  prod=ones(sn,1);
  for mu=1:n
    if mu~=nu
      prod=(xw0(:,1)-xw(mu,1)).*prod/(xw(nu,1)-xw(mu,1));
    end
  end
  prod=prod.^(2*s+1);
  for sig=0:2*s
    muhat(sig+1)=((xw0(:,1)'-xw(nu,1)).^sig.*prod')*xw0(:,2);
  end
  dtau=xw(:,1)-xw(nu,1); 
  for l=1:2*s
    u(l)=sum(1./dtau(find(dtau)).^l);
  end
  for k=1:2*s+1
    Ahat(k,k)=1;
  end
  for j=1:2*s
    for k=1:2*s+1-j
      Ahat(k,k+j)=-(2*s+1)*(u(1:j)*Ahat(1:j,j))/j;
    end
  end
  xw(nu,2:2*s+2)=(Ahat\muhat')';
  fac=1;
  for sig=0:2*s
    xw(nu,sig+2)=xw(nu,sig+2)/fac;
    fac=(sig+1)*fac;
  end
end   
