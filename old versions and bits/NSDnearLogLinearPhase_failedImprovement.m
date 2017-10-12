function [ X, W, R ] = NSDnearLogLinearPhase( N, k, a, b, m_g, c_g, SingPoint )
%steepest descent quadrature routine for linear phase functions
    %assumes g(x)=m_g*x + c_g
    
    %repeated computation
    osc=k*m_g;
    Rs=.15;  %radius inside of which log is considered to be nearly singular
    %p_max should be decreased in accordance with how far away the
    %singularity is.
    %grade towards singularity, then use Gauss-Laguerre along the rest of
    %the nearly-singular steepest descent paths
    if SingPoint < a
        [x0, w0]= genGaussLog(N,SingPoint,a,(a-SingPoint),'L' );
        %integral moves negatively along real axis, so negate the weights
        w0=-1*w0;
        %now do generalised Gauss Laguerre on (0,1]
        [xa, wa] = quad_gengauss_loglaguerre(N);
        %and warp these weights and nodes to [SingPoint, \infty]
        xa=xa+SingPoint;
        
        %get weights and nodes for second SD path
        [xb, wb] = GaussLaguerre(N, 0);
        xb = xb + b;
    elseif SingPoint > b
        [xa, wa] = GaussLaguerre(N, 0);
        x_a=x_a+a;
        
        [x0, w0]= genGaussLog(N,b,SingPoint,(SingPoint-b),'R' );
        %no need to negate the weights in this case
        [xb, wb] = quad_gengauss_loglaguerre(N);
        %and warp these weights and nodes to [SingPoint, \infty]
        xb=xb+SingPoint;
        
        %combine all weights and nodes for first SD path
        xb=[xb1; xb2;]; wb=[wb1; wb2];
    else
        error('Invalid singularity posish');
    end
    
    %weights
    W=(1i/osc)*exp(1i*k*c_g)*[exp(a*1i*osc)*wa;  -exp(b*1i*osc)*wb;];
    
    %nodes
    X1=1i*xa/osc;  X2=1i*xb/osc;
    X=[a+X1; b+X2;];
    %R=[sqrt(X1.^2 + (a-SingPoint)^2); sqrt(X2.^2 + (b-SingPoint)^2);];
    %R=[sqrt((xa/osc).^2 + (a-SingPoint)^2); sqrt((xb/osc).^2 + (b-SingPoint)^2);];
    if SingPoint < a
        R=X-SingPoint;
    else
        R=SingPoint-X;      
    end
    
end