function [ x, w ] = GaussHermiteTrunc( N,b )
%returns N weights and nodes for integral:
    % \int_0^b f(x)exp(-x^2) dx
    
    %there's a high chance that I've messed up with the k aligment here,
    %i.e. k+1 in the paper is equal to k here, but sometimes equal to k+1,
    %due to indexing from one, not zero
    
    N=N+1; %shift so can use vectors to store values
    
    gamma=NaN(N,1); alpha=NaN(N,1); beta=NaN(N,1);
    
    %have format of gamma(k+1):=\gamma_k, for k=0,1,...
    
    gamma(1) = 1/2*sqrt(pi)*erf(b);
    gamma(2) = 1/4*((2*exp(-2*b^2)*(exp(b^2) - 1)^2*(sqrt(pi) - 2*sqrt(pi)))/(pi*erf(b)) - 2*b*exp(-b^2) + sqrt(pi)*erf(b));
    indefInt1=@(x, k, alpha, beta)  (exp(-x.^2).*p_k(x, k, alpha, beta, b).*p_k(x, k-1, alpha, beta, b));
    indefInt2=@(x, k, alpha, beta) (exp(-x.^2).*(p_k(x, k, alpha, beta, b).^2));
    
    %get alpha and beta for base cases
    for k=2
        alpha(k) = .5/gamma(k) .* (indefInt2(b, k, alpha, beta)-indefInt2(0, k, alpha, beta));
        beta(k)  = -gamma(k)/gamma(k-1);
    end
    
    %now get higher gamma, which gives us beta and alpha
    for k=3:N
        gamma(k) = .5*(k-1)*gamma(k-1) - .5*(indefInt1(b, k, alpha, beta)-indefInt1(0, k, alpha, beta));
        alpha(k) = .5/gamma(k) .* (indefInt2(b, k, alpha, beta)-indefInt2(0, k, alpha, beta));
        beta(k)  = -gamma(k)/gamma(k-1);
    end
    
    % nodes are the roots of the function:
    p_N=@(t) p_k(t, N, alpha, beta, b);
    
    chebP_N=chebfun(@(t) p_N(t), [0 b]);
    
    x=roots(chebP_N);
    
    if length(x)~=(N-1)
        error('Something has gone wrong, there are too few roots for number of nodes required');
    end
    
    %weights given by
    w=gamma(N-1)./(dpdx_k(x, N, alpha, beta, b).*p_k(x, N-1, alpha, beta, b));
    
end