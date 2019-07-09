function Vout = forceDirichlet(x_1,x_2,Vin,maxDist,boundary)
%Forces Dirichlet condition in some neighbourhood of scatterer. A bit of a
%bodge to avoid computing the Hankel function with small arguments,
%and should (morally) only be used for making pretty plots.

    Vout = Vin;
    
    [X_1,X_2] = meshgrid(x_1,x_2);
    for m = 1:boundary.numComponents
        Nsamples = max(ceil(boundary.component(m).L/maxDist),100);
        t = linspace(0, boundary.component(m).L, Nsamples);
        Y = boundary.component(m).trace(t);
        Y1 = Y(1,:);
        Y2 = Y(2,:);
        for n = 1:Nsamples
            neighbInds = sqrt((X_1-Y1(n)).^2+(X_2-Y2(n)).^2)<maxDist;
            Vout(neighbInds) = 0;
        end
    end
    
end