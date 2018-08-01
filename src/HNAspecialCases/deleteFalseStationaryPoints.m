function [SPsOut, OrdersOut] = deleteFalseStationaryPoints(phaseDir,SPsIn, OrdersIn, thresh)
%sometimes the stationary points are false. This gets rid of them.
    if nargin == 3
        thresh = 1E-8;
    end
    SPsOut = [];
    OrdersOut = [];
    for n = 1:length(SPsIn)
        if abs(phaseDir(SPsIn(n))) < thresh
            SPsOut = [SPsOut SPsIn(n)];
            OrdersOut = [OrdersOut OrdersIn(n)];
        end
    end
end