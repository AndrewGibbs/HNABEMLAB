function [X,W] = abortPathSearch(a,b,ainf,binf,R,freq,g,Npts)
%if no path can be found, this file is called
    if ainf
        a=a*R;
    end
    if binf
        b=b*R;
    end
    [X,W_] = oscQuadExpensive(a,b,freq, Npts);
    W=W_.*exp(1i*freq*g(X));
end

