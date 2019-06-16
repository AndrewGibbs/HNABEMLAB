function P = p_kF( x, k, alpha, beta, b, p_power )
%
    if k==1
        P=ones(size(x));
    elseif k==2
        P=x-FreudMoment( b, p_power, 1 )/FreudMoment( b, p_power, 0 ); 
    else
        P=p_kF(x,k-1,alpha,beta,b,p_power).*(x-alpha(k-1)) - beta(k-1)*p_kF(x,k-2,alpha,beta,b,p_power);
    end

end