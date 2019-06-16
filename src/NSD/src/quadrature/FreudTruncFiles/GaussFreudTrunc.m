function [ x, w ] = GaussFreudTrunc( N, b, p)
%returns N weights and nodes for integral:
    % \int_0^b f(x)exp(-x^2) dx
    
    %there's a high chance that I've messed up with the k aligment here,
    %i.e. k+1 in the paper is equal to k here, but sometimes equal to k+1,
    %due to indexing from one, not zero
    
    N=N+1; %shift so can use vectors to store values
    
    gamma=NaN(N,1);
    
    %have format of gamma(k+1):=\gamma_k, for k=0,1,...
    
    moments=FreudMoment( b, p, 0:(2*N) );
    AlphaBeta=chebyshev(N+1,moments);
    % *probably got the indexing wrong here:
    alpha=[AlphaBeta(:,1)];
    beta=[AlphaBeta(:,2)];
    
    %now create the gamma vector (weighted norms of Gauss polynomials)
    gamma(1)=FreudMoment( b, p, 0 );
    for n=2:N
       gamma(n)=gamma(n-1)*beta(n);
    end
    
    % nodes are the roots of the function:
    p_N=@(t) p_kF(t, N, alpha, beta, b, p);
    
    chebP_N=chebfun(@(t) p_N(t), [0 b]);
    
    x=roots(chebP_N);
    
    if length(x)~=(N-1)
        error('Something has gone wrong, there are too few roots for number of nodes required');
    end
    
    %weights given by
    w=gamma(N-1)./(dpdx_kF(x, N, alpha, beta, b, p).*p_kF(x, N-1, alpha, beta, b, p));
    
end