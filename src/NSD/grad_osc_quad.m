function [X,W] = grad_osc_quad(a,b,ab_width,d,k,qppw,LR)
%quadrature rule designed for singular and/or oscillatory functions
    C1 = 2; %singularity parameter
    C2 = 2; %number of wavelengths in each subinterval
    delta = 0.15;
    max_layer = 15;
    num_layers = min(floor(C1*ab_width/d),max_layer);
    wavelength = 2*pi/k;
    [xG,wG]=gauleg(qppw);
    wG = (wG)/2;
    xG = ((xG)+1)/2;
    
    %first split interval into graded layers
    X = [];
    W = [];
    for n=0:num_layers
        %y = a + ab_width*delta^(num_layers - n);
        if n==0
            layer_width = ab_width*delta^num_layers;
            start_shift = 0;
        else
            layer_width = ab_width*delta^(num_layers-n)*(1-delta);
            start_shift = ab_width*delta^(num_layers-n+1);
        end
        num_splits = ceil(layer_width / (C2*wavelength));
        for m = 1:num_splits
            if strcmp(LR,'L')
                x = a + start_shift + (m-1)*layer_width/num_splits + xG*layer_width/num_splits;
            elseif strcmp(LR,'R')
                x = b - (start_shift + (m-1)*layer_width/num_splits + xG*layer_width/num_splits);
            end
            w = wG*layer_width/num_splits;
            X = [X; x];
            W = [W; w];
        end
    end
    if strcmp(LR,'R')
        X = flipud(X);
        W = flipud(W);
    end
end

