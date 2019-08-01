function vals = FarField(boundary, density, k, theta)
%produces matrix of domain values for density function
    theta = theta(:).';
    vals = zeros(length(theta),1);
    for n = 1:boundary.numComponents
        ffKernel = @(y1,y2,theta) exp(-1i*k*(cos(theta).*y1 + sin(theta).*y2));
        %first get values of density all around boundary
        t = density.edgeComponent(n).nodes;
        w = density.edgeComponent(n).weights;
        vt = density.edgeComponent(n).eval(t);
        Y = boundary.component(n).trace(t);
        Y1 = Y(:,1);
        Y2 = Y(:,2);
        vals = vals + (ffKernel(Y1,Y2,theta).*vt).'*w;
    end
end