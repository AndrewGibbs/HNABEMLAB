function vals = FarField(boundary, density, k, theta)
%produces matrix of domain values for density function

    ffKernel = @(y1,y2,theta) exp(-1i*k*(cos(theta).*y1 + sin(theta).*y2));
    %first get values of density all around boundary
    t = density.nodes;
    w = density.weights;
    vt = density.eval(t);
    Y = boundary.trace(t);
    Y1 = Y(:,1);
    Y2 = Y(:,2);
    vals = (ffKernel(Y1,Y2,theta).*vt).'*w;
end