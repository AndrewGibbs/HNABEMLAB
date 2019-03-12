function w = RiemannWeights(x,a,b)
%computes Riemann weights, as in examples in FNA2
    width = b-a;
    x_0 = x(end) - width;
    x_np1 = x(1) + width;
    xTended = [x_0; x; x_np1];
    w = .5*(xTended(3:end) - xTended(1:(end-2)));
end