function yn = constantPhase(a,b,G,Npts,thresh)

    [z, w] = quad_gauss(Npts, a, b);
    L2est = w.' * (G{2}(z).*conj(G{2}(z)));
    if L2est < thresh
        yn = true;
    else
        yn = false;
    end
end

