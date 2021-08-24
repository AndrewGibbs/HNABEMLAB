function domainPlot(Gamma,uinc,GOA,v_N,kwave,totalPixels)
    if nargin == 5
        totalPixels = 1000000;
    end
    V = Gamma.vertices;
    r = 1;
    xmin = min(V(:,1))- r; xmax = max(V(:,1)) + r;
    ymin = min(V(:,2)) - r; ymax = max(V(:,2)) + r;
    imageArea = (xmax-xmin)*(ymax-ymin);
    pixelRate = ceil(sqrt(totalPixels/imageArea));

    y = linspace(ymin,ymax,pixelRate*(ymax-ymin));
    x = linspace(xmin,xmax,pixelRate*(xmax-xmin));
    [X1, X2] = meshgrid(x,y);
    Sv = singleLayerDomain(Gamma, v_N, kwave, X1, X2);
    SPsi = singleLayerDomain(Gamma, GOA, kwave, X1, X2);
     u_N = uinc.eval(X1,X2) - Sv - SPsi;
%    u_N = Sv + SPsi;
    imagesc(x,y,real(u_N));
    set(gca,'YDir','normal');
    shading interp;
    hold on;
    axis equal;
    Gamma.draw;
end