function [gStationaryPointsOut, gSPordersOut] = pruneStationaryPoints(a, b, rectRad, gStationaryPointsIn, gSPordersIn, ainf, binf, settleRad)
%creates rectangle around neighbourhood of [a,b] and only keeps stationary
%points which are inside of it

    if isempty(gStationaryPointsIn)
        gStationaryPointsOut = [];
        gSPordersOut = [];
        return;
    end
    
     %This one should work better, although the '2' is quite arbitrary.
    if isempty(rectRad)
        if isinf(a+b)
            rectRad=settleRad;
        else
            rectRad=.75*(b-a);
        end
    end
    
    rect=[a-rectRad-rectRad*1i  b+rectRad-rectRad*1i  b+rectRad+rectRad*1i  a-rectRad+rectRad*1i]; 
    gStationaryPointsOut = [];
    gSPordersOut = [];
    for n = 1:length(gStationaryPointsIn)
        z = gStationaryPointsIn(n);
        if min(real(rect)) < real(z) && max(real(rect)) > real(z) && min(imag(rect)) < imag(z) && max(imag(rect)) > imag(z)
            gStationaryPointsOut = [gStationaryPointsOut z];
            gSPordersOut         = [gSPordersOut gSPordersIn(n)];
        end
    end
end

