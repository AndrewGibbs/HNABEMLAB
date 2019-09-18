function vals = FarField_lessSlow_stillSteady(boundary, density, k, theta)
%produces matrix of domain values for density function
    memConst = 15000^2; %based on experiments elsewhere
    
    %get width of matrix:
    N = 1;
    for n = 1:boundary.numComponents
        N = max(length(density.edgeComponent(n).nodes),N);
    end
    
    %get length of each one, given width and memory restrictions
    theta = theta(:).';
    thetaSegSize = ceil(memConst/(kwave*N));
    numThetaSegs = ceil(length(theta)/thetaSegSize);
    
    for m = 1:(numThetaSegs-1)
        thetaIndices{m} = (m-1)*thetaSegSize+1:(m*thetaSegSize);
        thetaSeg{m} = theta(thetaIndices{m});
        thetaLength(m) = length(thetaIndices{m});
    end
    
    % ADD THE FINAL BITs!
    thetaIndices{numThetaSegs} = (numThetaSegs*thetaSegSize):(length(theta));
    thetaSeg{numThetaSegs} = theta(thetaIndices{numThetaSegs});
    thetaLength(numThetaSegs) = length(thetaIndices{numThetaSegs});
    
    parfor m = 1:numThetaSegs
        valsSplit{m} = zeros(thetaLength,1);
        densityCpy = density;
        for n = 1:boundary.numComponents
            ffKernel = @(y1,y2,theta) exp(-1i*k*(cos(theta).*y1 + sin(theta).*y2));
            %first get values of density all around boundary
            t = densityCpy.edgeComponent(n).nodes;
            w = densityCpy.edgeComponent(n).weights;
            vt = densityCpy.edgeComponent(n).eval(t);
            Y = boundary.component(n).trace(t);
            Y1 = Y(:,1);
            Y2 = Y(:,2);
            valsSplit{m} = valsSplit{m} + (ffKernel(Y1,Y2,thetaSeg{m}).*vt).'*w;
        end
    end
    
    %now piece it back together
    vals = zeros(length(theta),1);
    for m = 1:numThetaSegs
        vals(thetaSeg{m}) = valsSplit{m};
    end
end