function [stationaryPoints,orders,badPoints] = symbolicStationaryPoints(x, fun, funSide, phase)
%looks at types of phase and determines if stationary points can be
%determined symbolically

    if nargin == 3
        tweaking = false;
    else
        tweaking =  true;
        thresh1 = 1E-8;
        thresh2 = .1;
    end

    Edge = fun.domain;
    %defaults, incase they can't be determined symbolically:
    stationaryPoints = NaN;
    orders = NaN;

    if isa(fun,'GeometricalOpticsFunction')
        if isa(fun.uinc,'planeWave')
            [stationaryPoints,orders] = HNApolygonSpecialCasePWGOA(x,fun.domain.side{funSide},fun.uinc.d);
            badPoints = nonAnalPoint(x,fun.domain.side{funSide});
        end
    elseif isa(fun,'baseFnHNA')
        [stationaryPoints,orders] = HNApolygonSpecialCase(x,Edge,fun.pm);
          badPoints = nonAnalPoint(x,Edge);
    else
        %(still plenty of room for expansion here)
        warning('Unable to determine stationary points for this type of phase');
    end
    
    
    
    %now remove rounding errors using a Newton method, if appropriate
    if tweaking
        stationaryPointsInit = stationaryPoints;
        stationaryPoints = [];
        for n = 1:length(stationaryPointsInit)
           if abs(phase{2}( stationaryPointsInit(n))) < thresh1
               stationaryPoints = [stationaryPoints stationaryPointsInit(n)];
           elseif  abs(phase{2}( stationaryPointsInit(n))) < thresh2 && orders(n)<2
               s = stationaryPointsInit(n);
               while abs(phase{2}(s)) >= thresh1
                    s = s - phase{2}(s)/phase{3}(s);
               end
               stationaryPoints = [stationaryPoints s];
           end
        end
    end
    
    %final check:
    for n = 1:length(stationaryPoints)
        if orders(n)>1
           for m = 2:orders(n)
               if abs(phase{m+1}(stationaryPoints(n))) > thresh1
                   orders(n) = orders(n) - 1;
                   break;
               end
           end
        end
    end
end