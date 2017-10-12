function [ X, W ] = GradedQuad( N,p_max,delta )
%1D graded quadrature routine
    [Z,w]=gauleg(N);
    Z=(Z+1)/2; w=w/2;
    X=delta^p_max*Z; W=delta^p_max*w;
    for p=fliplr(0:p_max-1)
        X=[X; delta^(p).*(delta+(1-delta)*Z)];
        W=[W; w*(1-delta)*delta^p];
    end
end