function [criticalPoints,pathPowers] = makeCriticalPoints(a,b,gStationaryPoints,spPowers,...
                                                             RectTol,gAnalytic, ainf, binf)
    
    %exclude cases where a and/or b are infinite points in complex plane
    if ainf
        a=[];
    end
    if binf
        b=[];
    end
    
    %now determine if a or b are at a stationary point
    if isempty(gStationaryPoints)
        criticalPoints=[a b];
        pathPowers=ones(1,~ainf+~binf);
    else
        criticalPoints=gStationaryPoints;
        pathPowers=spPowers;
        if ~isempty(a) && min(abs(gStationaryPoints-a))>RectTol
            criticalPoints=[a criticalPoints];
            pathPowers=[1 pathPowers];
        end
        if ~isempty(b) && min(abs(gStationaryPoints-b))>RectTol
            criticalPoints=[criticalPoints b];
            pathPowers=[pathPowers 1];
        end
    end
    if gAnalytic
        pathPowers=round(pathPowers);
    end
end