function P = dpdx_kF( x, k, alpha, beta, b, p )
%
    if k==1
        P=zeros(size(x));
    elseif k==2
        P=1; 
    else
        P=dpdx_kF(x,k-1,alpha,beta,b, p).*(x-alpha(k-1))  + p_kF(x,k-1,alpha,beta,b, p) - beta(k-1)*dpdx_kF(x,k-2,alpha,beta,b, p);
    end

end