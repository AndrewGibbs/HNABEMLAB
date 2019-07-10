function domainPlot(Gamma,uinc,GOA,v_N,kwave,totalPixels)
    if nargin == 5
        totalPixels = 1000000;
    end
    V = Gamma.vertices;
    xmin = min(V(1,:))-.3; xmax = max(V(1,:)) + .3;
    ymin = min(V(1,:)) -.2; ymax = max(V(2,:)) +.2;
    imageArea = (xmax-xmin)*(ymax-ymin);
    pixelRate = ceil(sqrt(totalPixels/imageArea));

    y = linspace(ymin,ymax,pixelRate*(ymax-ymin));
    x = linspace(xmin,xmax,pixelRate*(xmax-xmin));
    Sv = singleLayerDomain(Gamma, v_N, kwave, x, y);
    SPsi = singleLayerDomain(Gamma, GOA, kwave, x, y);
    [X1, X2] = meshgrid(x,y);
    u_N = uinc.eval(X1,X2) - Sv.' - SPsi.';
    imagesc(x,fliplr(y),flipud(real(u_N)));
    shading interp;
    set(gca,'YDir','normal');
    hold on;
    Gamma.draw;
end