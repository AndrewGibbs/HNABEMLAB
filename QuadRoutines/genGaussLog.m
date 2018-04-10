function [ X, W, R ] = genGaussLog( n,a,b,width,SingAt_LR )
%scaled generalised Gasusian quadrature for logarithmic singularity
    if b<a
        error('L<R');
    end
    if nargin<=3
        width=b-a;
    end
    if nargin<=4
        SingAt_LR='L';
    end
    
    %get standard weights and nodes over [0,1]
    [x,w]=quad_gengauss_log(n);
    
    if SingAt_LR=='L'
        X=a+width*x;
        W=width*w;
        R=width*x;
    elseif SingAt_LR=='R'
        X=b-width*flipud(x);
        W=flipud(width*w);
        R=width*flipud(x);
    else
        error('invalid entry for SingAt_LR');
    end

end