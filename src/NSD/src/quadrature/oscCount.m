function N = oscCount( a,b,freq,stationaryPoints,g )
%estimates the number of oscillations over [a,b], given stationary points,
%hence g is monotonic between them. If there are no stationary points, the
%user must input 'stationaryPoints=[]'

    s=[a stationaryPoints b];
    for j=1:(length(s)-1)
        N(j)=freq*abs(g(s(j+1))-g(s(j)));
    end

end

