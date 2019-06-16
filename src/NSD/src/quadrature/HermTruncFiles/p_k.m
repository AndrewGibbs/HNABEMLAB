function P = p_k( x, k, alpha, beta, b )
%
    if k==1
        P=ones(size(x));
    elseif k==2
        P=x-(1-exp(-b^2))/(sqrt(pi)*erf(b)); 
    else
        P=p_k(x,k-1,alpha,beta,b).*(x+alpha(k-1)) + beta(k-1)*p_k(x,k-2,alpha,beta,b);
    end

end