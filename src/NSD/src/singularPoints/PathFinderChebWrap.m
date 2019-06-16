% *** make sure ChebFun and PathFinder are added to path before running this
%andrew/Chebfun
function I = PathFinderChebWrap(a,b,freq,N,amp,Gphase,fSing, SP, ellipXtension)
%integrator for highly oscillatory functions which are analytic
%inside of a neighbourhood N([a,b]), but badly behaved outside of N([a,b])
    if nargin <=7
       SP = [];
       order = [];
    end
    if ~isempty(SP)
        order = 1;
    else
        order = [];
    end
    if nargin <= 8
        ellipXtension = 0;
    end
    minOscs = 30;
    Nchebs  = 20;
    N = max(ceil(Nchebs/2),N);
    %minOscs = 5;
    %new variable ellipXtension increases the interval of Chebyshev
    %approximation, to increase the size of the Bernstien ellipse, which
    %should improve the conditioning of the analytic extension in the
    %complex plane, such that IF 'amp' doesn't grow rapidly in
    %neighbourhood of [a,b], then neither will the analytic extension of
    %f|_{[a,b]}
    
    %now ChebFun it:
    AMP = chebfun(amp,[a-ellipXtension b+ellipXtension], Nchebs);
    chebPhase = chebfun( Gphase{1}, [a-ellipXtension b+ellipXtension], Nchebs);
    %DchebPhase = chebfun( Gphase{2}, [a b] );
    DchebPhase = diff(chebPhase);
    DDchebPhase = diff(DchebPhase);
    Gcheb = {@(x) chebPhase(x), @(x) DchebPhase(x), @(x) DDchebPhase(x)};%, @(x) DDchebPhase(x)};
 
    %now compute SD path, and evaluate integrals:
    [ z, w ] = PathFinder( a, b, freq, N, Gcheb, 'stationary points', SP, 'order', order,'fsingularities',fSing,'minoscs',minOscs);
    for path = [1 2]
        if max(abs(imag(z)))>0
            N = length(z)/2;
           if abs(AMP(w(path*N-1)*z(path*N-2))) < abs(w(path*N-1)*AMP(z(path*N-1))) && abs(w(path*N-1)*AMP(z(path*N-1))) > eps
               warning('Steepest descent path is not long enough, integrand still large at farthest node'); 
           end
        end
    end
    I = w.'*AMP(z);
end

