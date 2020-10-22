%getting singular values close to edges of domain - because of Hankel
%function. In particular a problem at low freqs, because range of
%singularity is much larger. I think.

function vals = singleLayerDomain(boundary, density, k, x1, x2)
%produces matrix of domain values for density function
    inshape = size(x1);
    dims_x1 = size(x1);
    dims_x2 = size(x2);
    if length(size(dims_x1)) == 2
        if dims_x1(1) ~= dims_x2(1) || dims_x1(2) ~= dims_x2(2)
            error('dimensions of grid do not match');
        end
        out_shape = [prod(dims_x1) 1];
    else
        out_shape = [dims_x1*dims_x2 1];
    end
    X1_ = reshape(x1,out_shape);
    X2_ = reshape(x2,out_shape);
    vals = zeros(inshape);
    Phi = @(x1,x2,y1,y2) 1i/4*besselh(0,1,k*sqrt((x1-y1).^2+(x2-y2).^2));
    %first get values of density all around boundary
    for n = 1:boundary.numComponents
        t = density.edgeComponent(n).nodes;
        w = density.edgeComponent(n).weights;
        vt = density.edgeComponent(n).eval(t);
        Y = boundary.component(n).trace(t);
        Y1 = Y(:,1);
        Y2 = Y(:,2);
        %now reshape these matrices into vectors
        Phi_R = Phi(X1_,X2_,Y1.',Y2.');
        vals_ = Phi_R*(vt.*w);
        vals = vals + reshape(vals_,inshape);
    end
end