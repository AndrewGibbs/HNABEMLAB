function intVal = PathFinderChebWrapGrad(sing,a,b,freq,N,amp,Gphase)
    %parameters for grading
    singPos = sing.position;
    sigma = 0.15;
    Nmax = 15;
    layerConst = 1;
    minOscs = 5;
    SP = []; order = [];
    
    %paramteres for ellipse extension:
    C1 = .5; % how close to singularity from interval
    C2 = 0; % how many interval lengths to extend interval by
    % (set to zero to zero to force no extension)
    
    %integrate the non-oscillatory bit first, if there is one
    x = findNonOscBit(Gphase{1},a,b,freq,minOscs);
    [ z, w ] = PathFinder( a, x, freq, N, Gphase, 'stationary points', SP, 'order', order,'fsingularities',sing,'minoscs',inf);
    J = w.'*amp(z);
    
    if b == x
        I = 0;
    else
        %decide how many layers of grading are appropz
        layers = min(Nmax,ceil(abs(layerConst*log(abs(singPos-x)/abs(b-x)))));
        %split [0,b] into subintervals:
        mesh(1) = x;
        mesh = [mesh  x+abs(b-x)*sigma.^((layers-1):-1:0)];

        for n = 1:layers
            %decide how much to extend the Bernstein ellipse by:
            r = min(abs(mesh(n)-singPos)*C1,abs(mesh(n+1)-mesh(n))*C2);
            I(n) = PathFinderChebWrap(mesh(n),mesh(n+1),freq,N,amp,Gphase,[],[],r);
            %... would be useful if widths could be passed too...
        end
    end
    intVal = J + sum(I);
end