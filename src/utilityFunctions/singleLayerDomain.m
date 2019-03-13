function vals = singleLayerDomain(boundary, density, k, x1, x2)
%produces matrix of domain values for density function

    Phi = @(x1,x2,y1,y2) 1i/4*besselh(0,1,k*sqrt((x1-y1).^2+(x2-y2).^2));
    %first get values of density all around boundary
    t = density.nodes;
    w = density.weights;
    vt = density.eval(t);
    Y = boundary.trace(t);
    Y1 = Y(:,1);
    Y2 = Y(:,2);
    [X1,X2] = meshgrid(x1,x2);
    %now reshape these matrices into vectors
    [m,n]= size(X1);
    X1_ = reshape(X1,[m*n 1]);
    X2_ = reshape(X2,[m*n 1]);
    Phi_R = Phi(X1_,X2_,Y1.',Y2.');
    vals_ = Phi_R*(vt.*w);
    vals = reshape(vals_,[m n]);
end