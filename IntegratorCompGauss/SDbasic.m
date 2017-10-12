function [ X, W ] = SDbasic( N, k,a,b,m_g,c_g )
%steepest descent quadrature routine for linear phase functions
    %assumes g(x)=m_g*x + c_g
    %c_g=0;
    %repeated computation
    osc=k*m_g;

    %use Gauss Laguerre quadrature along steepest descent paths
    [x, w] = GaussLaguerre(N, 0);
    
    %weights
    W=(1i/osc)*exp(1i*k*c_g)*[exp(a*1i*osc)*w;  -exp(b*1i*osc)*w;];
    
    %for a bit of extra efficiency
    X_=1i*x/osc;
    %nodes
    X=[a+X_; b+X_;];
    
end