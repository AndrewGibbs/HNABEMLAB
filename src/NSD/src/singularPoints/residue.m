function [fDiff, z, W] = residue( f,z0,radius, N )
%computes the residue of a function with a singularity at z0
    if nargin<=5
        N=15;
    end

    [theta, w]=quad_gauss(N, 0, 2*pi);
    z=z0+radius*exp(1i*theta);
    
    %I thought the trapezium rule would work well, but it didn't. So using
    %Gaussian quadrature insstead.
    fDiff=w.'*(1i*radius*exp(1i*theta).*f(z));
    
    %create 'one' to divide by at the end
    %one=w.'*(1i*radius*exp(1i*theta).*(1./((z-centre))));
    one=1;
    
    %now transform output into same shape as input
    fDiff=fDiff/one;
    
    %output weights too, so this can all be combined with other weights and
    %nodes:
    W=w.*(1i*radius*exp(1i*theta))/one;
    
end