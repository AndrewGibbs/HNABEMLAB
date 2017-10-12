function [ X, W, R ] = NSDnearLogLinearPhase( N, k, a, b, m_g, c_g, SingPoint, p_max, GradDelta )
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
        pathSplit=sqrt(Rs^2-(SingPoint-a)^2);
        [xa1, wa1]= GradedQuad( N,p_max,GradDelta );
        %now scale the graded weights and nodes from [0 1] to [0 pathSplit]
        xa1=xa1*pathSplit; wa1=wa1*pathSplit.*exp(-xa1);
        %now do Gauss Laguerre on [pathSplit \infty]
        [xa2, wa2] = GaussLaguerre(N, 0);
        xa2=xa2+pathSplit; wa2=wa2*exp(-pathSplit);
        %combine all weights and nodes for first SD path
        xa=[xa1; xa2;]; wa=[wa1; wa2];
        
        %get weights and nodes for second SD path
        [xb, wb] = GaussLaguerre(N, 0);
    elseif SingPoint > b
        pathSplit=sqrt(Rs^2-(b-SingPoint)^2);
        [xa, wa] = GaussLaguerre(N, 0);
        
        [xb1, wb1]= GradedQuad( N,p_max,GradDelta );
        %now scale the graded weights and nodes from [0 1] to [0 pathSplit]
        xb1=xb1*pathSplit; wb1=wb1*pathSplit.*exp(-xb1);
        %now do Gauss Laguerre on [a \infty]
        [xb2, wb2] = GaussLaguerre(N, 0);
        xb2=xb2+pathSplit; wb2=wb2*exp(-pathSplit);
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