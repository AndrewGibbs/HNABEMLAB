function P = dpdx_k( x, k, alpha, beta, b )
%
    if k==1
        P=zeros(size(x));
    elseif k==2
        P=1; 
    else
        P=dpdx_k(x,k-1,alpha,beta,b).*(x+alpha(k-1))  + p_k(x,k-1,alpha,beta,b) + beta(k-1)*dpdx_k(x,k-2,alpha,beta,b);
    end

end