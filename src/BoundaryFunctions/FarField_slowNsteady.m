function vals = FarField_slowNsteady(boundary, density, k, theta)
%produces matrix of domain values for density function.
%more effective when there are a large number of observation points, avoids
%matrix computations so is memory efficient
    theta = theta(:).';
    vals = zeros(length(theta),1);
    ffKernel = @(y1,y2,theta) exp(-1i*k*(cos(theta).*y1 + sin(theta).*y2));
    parfor obsCount = 1:length(theta)
        for n = 1:boundary.numComponents
            %first get values of density all around boundary
            t = density.edgeComponent(n).nodes;
            w = density.edgeComponent(n).weights;
            vt = density.edgeComponent(n).eval(t);
            Y = boundary.component(n).trace(t);
            Y1 = Y(:,1);
            Y2 = Y(:,2);
            vals(obsCount) = sum(ffKernel(Y1,Y2,theta(obsCount)).*vt.*w);
        end
    end
end