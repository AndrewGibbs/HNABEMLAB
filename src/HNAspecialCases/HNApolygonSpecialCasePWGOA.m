function [stationaryPoints,orders] = HNApolygonSpecialCasePWGOA(x,side,d)
%computes (theotetically) the stationary point of the phase which is the
%analytic continuation of g(z) = |x-y(z)|+d.y(z)

    %get parametric parameters
    ydot = side.dSv;
    P = side.P1;
    
    %use incredibly long symbolic result to determine stationary points
    stationaryPointsInit(1) = (P(1)*ydot(1)^3 + P(2)*ydot(2)^3 - x(1)*ydot(1)^3 - x(2)*ydot(2)^3 + P(1)*ydot(1)*ydot(2)^2 + P(2)*ydot(1)^2*ydot(2) - x(1)*ydot(1)*ydot(2)^2 - x(2)*ydot(1)^2*ydot(2) - P(1)*d(1)^2*ydot(1)^3 - P(2)*d(2)^2*ydot(2)^3 + d(1)^2*x(1)*ydot(1)^3 + d(2)^2*x(2)*ydot(2)^3 + P(2)*d(1)*ydot(1)^2*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - P(1)*d(2)*ydot(2)^2*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - d(1)*x(2)*ydot(1)^2*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + d(2)*x(1)*ydot(2)^2*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - P(1)*d(2)^2*ydot(1)*ydot(2)^2 - P(2)*d(1)^2*ydot(1)^2*ydot(2) + d(1)^2*x(2)*ydot(1)^2*ydot(2) + d(2)^2*x(1)*ydot(1)*ydot(2)^2 - P(1)*d(1)*ydot(1)*ydot(2)*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + P(2)*d(2)*ydot(1)*ydot(2)*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + d(1)*x(1)*ydot(1)*ydot(2)*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - d(2)*x(2)*ydot(1)*ydot(2)*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - 2*P(1)*d(1)*d(2)*ydot(1)^2*ydot(2) - 2*P(2)*d(1)*d(2)*ydot(1)*ydot(2)^2 + 2*d(1)*d(2)*x(1)*ydot(1)^2*ydot(2) + 2*d(1)*d(2)*x(2)*ydot(1)*ydot(2)^2)/(d(1)^2*ydot(1)^4 + d(1)^2*ydot(1)^2*ydot(2)^2 + 2*d(1)*d(2)*ydot(1)^3*ydot(2) + 2*d(1)*d(2)*ydot(1)*ydot(2)^3 + d(2)^2*ydot(1)^2*ydot(2)^2 + d(2)^2*ydot(2)^4 - ydot(1)^4 - 2*ydot(1)^2*ydot(2)^2 - ydot(2)^4);
    stationaryPointsInit(2) = (P(1)*ydot(1)^3 + P(2)*ydot(2)^3 - x(1)*ydot(1)^3 - x(2)*ydot(2)^3 + P(1)*ydot(1)*ydot(2)^2 + P(2)*ydot(1)^2*ydot(2) - x(1)*ydot(1)*ydot(2)^2 - x(2)*ydot(1)^2*ydot(2) - P(1)*d(1)^2*ydot(1)^3 - P(2)*d(2)^2*ydot(2)^3 + d(1)^2*x(1)*ydot(1)^3 + d(2)^2*x(2)*ydot(2)^3 - P(2)*d(1)*ydot(1)^2*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + P(1)*d(2)*ydot(2)^2*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + d(1)*x(2)*ydot(1)^2*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - d(2)*x(1)*ydot(2)^2*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - P(1)*d(2)^2*ydot(1)*ydot(2)^2 - P(2)*d(1)^2*ydot(1)^2*ydot(2) + d(1)^2*x(2)*ydot(1)^2*ydot(2) + d(2)^2*x(1)*ydot(1)*ydot(2)^2 + P(1)*d(1)*ydot(1)*ydot(2)*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - P(2)*d(2)*ydot(1)*ydot(2)*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - d(1)*x(1)*ydot(1)*ydot(2)*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + d(2)*x(2)*ydot(1)*ydot(2)*(- d(1)^2*ydot(1)^2 - 2*d(1)*d(2)*ydot(1)*ydot(2) - d(2)^2*ydot(2)^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - 2*P(1)*d(1)*d(2)*ydot(1)^2*ydot(2) - 2*P(2)*d(1)*d(2)*ydot(1)*ydot(2)^2 + 2*d(1)*d(2)*x(1)*ydot(1)^2*ydot(2) + 2*d(1)*d(2)*x(2)*ydot(1)*ydot(2)^2)/(d(1)^2*ydot(1)^4 + d(1)^2*ydot(1)^2*ydot(2)^2 + 2*d(1)*d(2)*ydot(1)^3*ydot(2) + 2*d(1)*d(2)*ydot(1)*ydot(2)^3 + d(2)^2*ydot(1)^2*ydot(2)^2 + d(2)^2*ydot(2)^4 - ydot(1)^4 - 2*ydot(1)^2*ydot(2)^2 - ydot(2)^4);
    
    %check for degenerates, and group nearby stationary points if there are
    %two (I guess there usually should be two)
    [stationaryPoints,orders]  = sortStationaryPoints(stationaryPointsInit,1E-8);
    
end