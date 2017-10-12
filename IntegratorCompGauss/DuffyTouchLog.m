%  Calculate Quadraturepoints for integrals of the form:
%
%  \int_{\tau_1}\int_{\tau_2} f(x,y) dxdy
%
%  where the function f has a singularity at [0,1]
%
%  INPUT:
%    n_gl  : the number of gauss legendre points
%    n_cgl : number for the composite gauss legendre quadrature
%               number of quadrature points is (n_cgl-1)*(n_cgl+2)/2+1
%
%  OUTPUT:
%  [X,W]   : Quadrature points on [0,1]^2 and quadrature weights
%
function [x,y,r,w, bmx, ymb]=DuffyTouchLog(a,b,c,angle,qppw,Grad_delta,near,Nlayers,ab_width,bc_width)
%this will need to adjust when close to the corners of the scatterer.

if nargin==8
    ab_width=b-a; bc_width=c-b;
end
if min(ab_width/bc_width,bc_width/ab_width)<near;
    if ab_width<bc_width
        delta=ab_width/bc_width;
    else
        delta=bc_width/ab_width;
    end
     near_max_p=min(ceil(log_base(delta,Grad_delta))+1,Nlayers);
     [tau,w_tau] = GradedQuad( qppw, near_max_p, Grad_delta );
else
    [tau,w_tau]=gauleg(qppw+1); tau=tau/2+1/2; w_tau=w_tau/2;
end

    %get the points
    [z,w_z]=quad_gengauss_log(qppw);
    [Z,TAU] = rearrange( z,tau );
    [W_Z,W_TAU] = rearrange( w_z,w_tau );
    
    %do the transform
    WZ=W_Z.*W_TAU.*Z; WZ=[WZ.' WZ.'].';

    %new arguments for the kernel
    X=[Z.' (Z.*TAU).'].';
    Y=[(Z.*TAU).' Z.'].';
    w=WZ*(ab_width*bc_width);
    
    %compute r:=|x-y| minimising rounding error, accounting for any angle
    %between mesh elements... 
    if angle~=pi
        r = sqrt((ab_width*X).^2 + (bc_width*Y).^2 ...
            -2*cos(angle)*(ab_width*X).*(bc_width*Y));
    else
        r=ab_width*X + bc_width*Y;
    end
    
    bmx=ab_width*X; ymb=bc_width*Y;
    x=b-bmx;        y=b+ymb;
end