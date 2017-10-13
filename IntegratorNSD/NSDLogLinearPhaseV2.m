function [ X, W ] = NSDLogLinearPhaseV2( N, k,a,b,m_g,c_g,SingPoint, p_max,GradDelta  )
%steepest descent quadrature routine for linear phase functions
    %assumes g(x)=m_g*x + c_g
    %repeated computation
    if nargin<=6
        SingPoint=inf;
    end
    osc=k*m_g;
    Rs=.15;
    %use Gauss Laguerre quadrature along steepest descent paths
    %first determine path from point 'a'
    if a==SingPoint
        [xa, wa]=quad_gengauss_loglaguerre(N);
    elseif abs(a-SingPoint)<=Rs %nearly singular
        pathSplit=sqrt(Rs^2-(SingPoint-a)^2);
        [xa1, wa1]= GradedQuad( N,p_max,GradDelta );
        %now scale the graded weights and nodes from [0 1] to [0 pathSplit]
        xa1=xa1*pathSplit; wa1=wa1*pathSplit.*exp(-xa1);
        %now do Gauss Laguerre on [pathSplit \infty]
        [xa2, wa2] = GaussLaguerre(N, 0);
        xa2=xa2+pathSplit; wa2=wa2*exp(-pathSplit);
        %combine all weights and nodes for first SD path
        xa=[xa1; xa2;]; wa=[wa1; wa2];
    else    %not even nearly singular
        [xa, wa] = GaussLaguerre(N, 0);
    end
    if b==SingPoint
        [xb, wb]=quad_gengauss_loglaguerre(N);
    elseif abs(b-SingPoint)<=Rs %nearly singular
        pathSplit=sqrt(Rs^2-(b-SingPoint)^2);
        [xb1, wb1]= GradedQuad( N,p_max,GradDelta );
        %now scale the graded weights and nodes from [0 1] to [0 pathSplit]
        xb1=xb1*pathSplit; wb1=wb1*pathSplit.*exp(-xb1);
        %now do Gauss Laguerre on [a \infty]
        [xb2, wb2] = GaussLaguerre(N, 0);
        xb2=xb2+pathSplit; wb2=wb2*exp(-pathSplit);
        %combine all weights and nodes for first SD path
        xb=[xb1; xb2;]; wb=[wb1; wb2];
    else %not even nearly singular
        [xb, wb] = GaussLaguerre(N, 0);
    end
    
    %weights
    W=(1i/osc)*exp(1i*k*c_g)*[exp(a*1i*osc)*wa;  -exp(b*1i*osc)*wb;];
    
    %nodes
    X=[a+1i*xa/osc; b+1i*xb/osc;];
    
end