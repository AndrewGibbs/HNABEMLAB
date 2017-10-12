function [ X, W, R ] = NSDLogLinearPhase( N, k,a,b,m_g,c_g,SingAt_LR )
%steepest descent quadrature routine for linear phase functions
    %assumes g(x)=m_g*x + c_g
    %repeated computation
    osc=k*m_g;

    %use Gauss Laguerre quadrature along steepest descent paths
    if strcmp(SingAt_LR,'L')
        [xa, wa]=quad_gengauss_loglaguerre(N);
        [xb, wb] = GaussLaguerre(N, 0);
    elseif strcmp(SingAt_LR,'R')
        [xa, wa] = GaussLaguerre(N, 0);
        [xb, wb]=quad_gengauss_loglaguerre(N);
    else
        error('Invalid singularity posish');
    end
    
    %weights
    W=(1i/osc)*exp(1i*k*c_g)*[exp(a*1i*osc)*wa;  -exp(b*1i*osc)*wb;];
    
    %nodes
    X=[a+1i*xa/osc; b+1i*xb/osc;];
%     if strcmp(SingAt_LR,'L')
%         R= [1i*xa/osc; abs(a-(b+1i*xb/osc));];
%     elseif strcmp(SingAt_LR,'R')
%         R= [abs(b-(a+1i*xa/osc)); 1i*xb/osc;];
%     end
    if strcmp(SingAt_LR,'L')
        R= X-a;
    elseif strcmp(SingAt_LR,'R')
        R= b-X;
    end
    
end