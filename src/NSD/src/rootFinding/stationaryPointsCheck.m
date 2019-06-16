function [gStationaryPointsOut, ordersOut] = stationaryPointsCheck(G,gStationaryPoints,orders)
    dgdx = G{2};
    thresh = 1E-8;
    gStationaryPointsOut = [];
    ordersOut = [];
    for n = 1:length(gStationaryPoints)
        if abs(dgdx(gStationaryPoints(n)))>thresh
            warning(sprintf('%f+%fi is not a stationary point of phase provided',real(gStationaryPoints(n)),imag(gStationaryPoints(n))));
        else
            gStationaryPointsOut = [gStationaryPointsOut gStationaryPoints(n)];
            ordersOut = [ ordersOut orders(n)];
        end
    end
end