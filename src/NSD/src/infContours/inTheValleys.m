function valSum = inTheValleys(a, g, R, Npts, freq)
%inverts region outside circle to inside of circle, and checks how well it
%approximates the true SD path to infinity
    if nargin<=3
        Npts=15;
    end
    if nargin<=4
        freq=1; %extra freq info can only speed up process
    end
    
    %order of monomial
    order=length(a)-1;
    
    %now determine valleys at infnity
    val = exp(2i*pi*(1/4 + (1:order))/order);
    
    %get some quadrature points, use them a few times
    [z,w] = quad_gauss(Npts, 0, 1);
    
    valSum = 0;
    for v=val
       valSum = valSum + abs ( (-1/R) * w.' * (exp(-freq*imag(g(v*R./z))).*(z.^2)) );
    end
end

