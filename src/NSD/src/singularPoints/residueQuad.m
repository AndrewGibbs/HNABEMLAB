function [z, W] = residueQuad( z0,radius, N )
%computes the residue of a function with a singularity at z0
    if nargin<=5
        N=15;
    end

    [theta, w]=quad_gauss(N, 0, 2*pi);
    z=z0+radius*exp(1i*theta);
    
    W=w.*(1i*radius*exp(1i*theta));
    
end