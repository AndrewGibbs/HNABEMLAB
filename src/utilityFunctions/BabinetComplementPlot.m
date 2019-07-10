function BabinetComplementPlot(Gamma,uinc,GOA,v_N,kwave,totalPixels)
    if nargin == 5
        totalPixels = 1000000;
    end
     V = Gamma.vertices;
    xmin = min(V(1,:))-.3; xmax = max(V(1,:)) + .3;
    ymin = min(V(1,:)) -.2; ymax = max(V(2,:)) +.2;
    imageArea = (xmax-xmin)*(ymax-ymin);
    pixelRate = ceil(sqrt(totalPixels/imageArea));
    pixelRate = ceil(sqrt(totalPixels/imageArea));
    
    yU = linspace(0,ymax,pixelRate*(ymax));
    yL = linspace(ymin,0,pixelRate*(-ymin));
    yL = yL(1:(end-1));
    x = linspace(xmin,xmax,pixelRate*(xmax-xmin));
    SvU = singleLayerDomain(Gamma, v_N, kwave, x, yU);
    SPsiU = singleLayerDomain(Gamma, GOA, kwave, x, yU);
    SvL = singleLayerDomain(Gamma, v_N, kwave, x, yL);
    SPsiL = singleLayerDomain(Gamma, GOA, kwave, x, yL);

    uincR = uinc.getReflection(Gamma);
    [X1U, X2U] = meshgrid(x,yU);
    u_NU = uinc.eval(X1U,X2U) - SvU.' - SPsiU.';

    DirZone = .01/kwave;
    B_U = uincR.eval(X1U,X2U) + forceDirichlet(x,yU,u_NU,DirZone,Gamma);

    [X1L, X2L] = meshgrid(x,yL);
    u_NL = uinc.eval(X1L,X2L) - SvL.' - SPsiL.';
    B_L = uinc.eval(X1L,X2L) - forceDirichlet(x,yL,u_NL,DirZone,Gamma);

    babz = [B_L; B_U;];
    imagesc(x,fliplr([yL yU]),flipud(real(babz)));

    shading interp;
    set(gca,'YDir','normal');
    hold on;
    Gamma.drawCompliment(xmin,xmax);

end

