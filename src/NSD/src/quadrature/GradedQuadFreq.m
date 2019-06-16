function [ X, W ] = GradedQuadFreq( N,p_max,delta, oscs )
%1D graded quadrature routine

    oscScaler=ceil(oscs*(1-delta));
    [Z,w]=gauleg(N*oscScaler);
    Z=(Z+1)/2; w=w/2;
    X=delta^p_max*Z; W=delta^p_max*w;
        
    for p=fliplr(0:p_max-1)
        oscScaler=ceil(oscs*(1-delta)^(p+1));
        [Z,w]=gauleg(N*oscScaler);
        Z=(Z+1)/2; w=w/2;
        X=[X; delta^(p).*(delta+(1-delta)*Z)];
        W=[W; w*(1-delta)*delta^p];
    end
end