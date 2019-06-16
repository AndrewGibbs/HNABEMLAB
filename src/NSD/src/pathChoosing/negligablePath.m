function tf = negligablePath( p1,p2, g,freq, thresh, Npts )
%compute line integral between two points, to decide if they can be
%considered connected in the context of the SD adjacency matrix
    %basicPath=linspace(p1,p2,Npts);
    [z,w] = quad_gauss(Npts, p1, p2);

    if abs(w.'*exp(-freq*imag(g(z)))) < thresh
        tf=true;
    else
        tf=false;
    end
    
end