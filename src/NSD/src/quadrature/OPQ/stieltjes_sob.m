% STIELTJES_SOB Stieltjes algorithm for Sobolev orthogonal polynomials.
%
%    The call [B,normsq]=STIELTJES_SOB(n,s,nd,xw,a0,same) generates
%    the nxn matrix B of the recurrence coefficients of the monic
%    Sobolev orthogonal polynomials, with beta_j^k, 0<=j<=k, 
%    k=0,1,...,n-1, occupying the position (j+1,k+1) of the matrix B
%    and all remaining elements of B being zero. The routine also
%    returns the n-vector normsq of the squared norms of the
%    Sobolev orthogonal polynomials.
%
%    The routine uses an inner-product representation of the
%    coefficients beta_j^k and evaluates the inner products exactly
%    by Gaussian quadrature rules. The number of points n_sigma in
%    these rules may be different for the s+1 individual measures
%    making up the Sobolev inner product. The input variable
%    nd collects the numbers n_sigma in a vector of dimension s+1.
%    The input variable xw is an array of dimension Nx(2s+2), where
%    N=max n_sigma. It stores the nodes of the quadrature rules
%    in the first half, and the weights in the second half of the
%    array. The input parameter a0 is the first alpha-coefficient
%    of the "ground measure" of the Sobolev inner product, and the
%    variable same has to be set equal to 1 if all quadrature rules
%    have the same set of nodes, and equal to 0 otherwise.
%
%    The procedure is an extended version of the Stieltjes algorithm
%    described in Section 4 of the paper cited below.
%
%    REFERENCE: W. Gautschi and M. Zhang, "Computing orthogonal
%    polynomials in Sobolev spaces", Numer. Math. 71 (1995), 259-183.
%
function [B,normsq]=stieltjes_sob(n,s,nd,xw,a0,same)
if n<1, error('n out of range'), end
if s<0, error('s out of range'), end
if s>=n, s=n-1; end
if any(nd<1), error('nd has nonpositive components'), end
md=max(nd); if size(xw,1)<md, error('columns of xw too short'), end
B=zeros(n); 
%
% Initialization
%
p=zeros(n,s+1,md); p(1,1,:)=1; xp=zeros(n,s+1,s+1,md);
for is=1:s+1
  for js=1:s+1
    if(same&is==js | ~same)
      for mu=1:nd(js)
        if is==1
          xp(1,is,js,mu)=xw(mu,js);
        elseif is==2
          xp(1,is,js,mu)=1;
        else
          xp(1,is,js,mu)=0;
        end
      end
    end
  end
end
B(1,1)=a0;
if n==1, return, end
%
% Continuation
%
for k=2:n
  for j=1:k
    for is=1:s+1
      for nu=1:nd(is)
        sm=xp(k-1,is,is,nu);
        for l=1:k-1
          sm=sm-B(l,k-1)*p(k-l,is,nu);
        end
        p(k,is,nu)=sm;
      end
    end
    Bnum=0;
    Bden=0;
    for is=1:s+1
      for js=1:s+1
        if(same&is==js | ~same)
          for mu=1:nd(js)
            xsm=xw(mu,js)*xp(k-1,is,js,mu);
            if is>1
              if same
                xsm=xsm+(is-1)*xp(k-1,is-1,is-1,mu);
              else
                xsm=xsm+(is-1)*xp(k-1,is-1,js,mu);
              end
            end
            for l=1:k-1
              xsm=xsm-B(l,k-1)*xp(k-l,is,js,mu);
            end
            xp(k,is,js,mu)=xsm;
            if is==js
              Bnum=Bnum+xw(mu,is+s+1)*xp(k,is,js,mu)*p(k-j+1,is,mu);
              Bden=Bden+xw(mu,is+s+1)*p(k-j+1,is,mu)^2;
            end
          end
        end
      end
    end
    B(j,k)=Bnum/Bden;
    if k==n, normsq(n-j+1)=Bden; end
  end
end
