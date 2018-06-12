function [stationaryPoints,orders] = HNApolygonSpecialCase(x,side,pm)
%computes (theotetically) the stationary point of the phase which is the
%analytic continuation of g(z) = |x-y(z)|+z
    
    %get parametric parameters
    ydot = side.dSv;
    P = side.P1;
    thresh = 1E-8;
    
    %use incredibly long symbolic result to determine stationary points
%     stationaryPointsInit(1) = -(P(1)*ydot(1)^3 + P(2)*ydot(2)^3 - x(1)*ydot(1)^3 - x(2)*ydot(2)^3 + P(1)*ydot(1)*ydot(2)^2 + P(2)*ydot(1)^2*ydot(2) + pm^2*x(1)*ydot(1) + pm^2*x(2)*ydot(2) - x(1)*ydot(1)*ydot(2)^2 - x(2)*ydot(1)^2*ydot(2) - P(1)*pm^2*ydot(1) - P(2)*pm^2*ydot(2) + P(1)*pm*ydot(2)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - P(2)*pm*ydot(1)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - pm*x(1)*ydot(2)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + pm*x(2)*ydot(1)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2))/(- pm^2*ydot(1)^2 - pm^2*ydot(2)^2 + ydot(1)^4 + 2*ydot(1)^2*ydot(2)^2 + ydot(2)^4);
%     stationaryPointsInit(2) = -(P(1)*ydot(1)^3 + P(2)*ydot(2)^3 - x(1)*ydot(1)^3 - x(2)*ydot(2)^3 + P(1)*ydot(1)*ydot(2)^2 + P(2)*ydot(1)^2*ydot(2) + pm^2*x(1)*ydot(1) + pm^2*x(2)*ydot(2) - x(1)*ydot(1)*ydot(2)^2 - x(2)*ydot(1)^2*ydot(2) - P(1)*pm^2*ydot(1) - P(2)*pm^2*ydot(2) - P(1)*pm*ydot(2)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + P(2)*pm*ydot(1)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + pm*x(1)*ydot(2)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - pm*x(2)*ydot(1)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2))/(- pm^2*ydot(1)^2 - pm^2*ydot(2)^2 + ydot(1)^4 + 2*ydot(1)^2*ydot(2)^2 + ydot(2)^4);
%     stationaryPointsInit(1) = (P(1)*ydot(1) + P(2)*ydot(2) - x(1)*ydot(1) - x(2)*ydot(2) - P(1)*ydot(1)^3 - P(2)*ydot(2)^3 + x(1)*ydot(1)^3 + x(2)*ydot(2)^3 - P(1)*ydot(1)*ydot(2)^2 - P(2)*ydot(1)^2*ydot(2) + x(1)*ydot(1)*ydot(2)^2 + x(2)*ydot(1)^2*ydot(2) + P(1)*ydot(2)*(ydot(1)^2 + ydot(2)^2 - 1)^(1/2) - P(2)*ydot(1)*(ydot(1)^2 + ydot(2)^2 - 1)^(1/2) - x(1)*ydot(2)*(ydot(1)^2 + ydot(2)^2 - 1)^(1/2) + x(2)*ydot(1)*(ydot(1)^2 + ydot(2)^2 - 1)^(1/2))/(ydot(1)^4 + 2*ydot(1)^2*ydot(2)^2 - ydot(1)^2 + ydot(2)^4 - ydot(2)^2);
%     stationaryPointsInit(2) = (P(1)*ydot(1) + P(2)*ydot(2) - x(1)*ydot(1) - x(2)*ydot(2) - P(1)*ydot(1)^3 - P(2)*ydot(2)^3 + x(1)*ydot(1)^3 + x(2)*ydot(2)^3 - P(1)*ydot(1)*ydot(2)^2 - P(2)*ydot(1)^2*ydot(2) + x(1)*ydot(1)*ydot(2)^2 + x(2)*ydot(1)^2*ydot(2) - P(1)*ydot(2)*(ydot(1)^2 + ydot(2)^2 - 1)^(1/2) + P(2)*ydot(1)*(ydot(1)^2 + ydot(2)^2 - 1)^(1/2) + x(1)*ydot(2)*(ydot(1)^2 + ydot(2)^2 - 1)^(1/2) - x(2)*ydot(1)*(ydot(1)^2 + ydot(2)^2 - 1)^(1/2))/(ydot(1)^4 + 2*ydot(1)^2*ydot(2)^2 - ydot(1)^2 + ydot(2)^4 - ydot(2)^2);
%     stationaryPointsInit(1) = -(P(1)*ydot(1)^3 + P(2)*ydot(2)^3 - x1*ydot(1)^3 - x2*ydot(2)^3 + P(1)*ydot(1)*ydot(2)^2 + P(2)*ydot(1)^2*ydot(2) + pm^2*x1*ydot(1) + pm^2*x2*ydot(2) - x1*ydot(1)*ydot(2)^2 - x2*ydot(1)^2*ydot(2) - P(1)*pm^2*ydot(1) - P(2)*pm^2*ydot(2) + P(1)*pm*ydot(2)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - P(2)*pm*ydot(1)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - pm*x1*ydot(2)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + pm*x2*ydot(1)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2))/(- pm^2*ydot(1)^2 - pm^2*ydot(2)^2 + ydot(1)^4 + 2*ydot(1)^2*ydot(2)^2 + ydot(2)^4)
%     stationaryPointsInit(2) = -(P(1)*ydot(1)^3 + P(2)*ydot(2)^3 - x(1)*ydot(1)^3 - x(2)*ydot(2)^3 + P(1)*ydot(1)*ydot(2)^2 + P(2)*ydot(1)^2*ydot(2) + pm^2*x(1)*ydot(1) + pm^2*x(2)*ydot(2) - x(1)*ydot(1)*ydot(2)^2 - x(2)*ydot(1)^2*ydot(2) - P(1)*pm^2*ydot(1) - P(2)*pm^2*ydot(2) - P(1)*pm*ydot(2)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + P(2)*pm*ydot(1)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) + pm*x(1)*ydot(2)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2) - pm*x(2)*ydot(1)*(- pm^2 + ydot(1)^2 + ydot(2)^2)^(1/2))/(- pm^2*ydot(1)^2 - pm^2*ydot(2)^2 + ydot(1)^4 + 2*ydot(1)^2*ydot(2)^2 + ydot(2)^4)
 
    %stationaryPoints = (P(1)^2*pm^2 - P(1)^2*ydot(1)^2 - 2*P(1)*P(2)*ydot(1)*ydot(2) - 2*P(1)*pm^2*x(1) + 2*P(1)*x(1)*ydot(1)^2 + 2*P(1)*x(2)*ydot(1)*ydot(2) + P(2)^2*pm^2 - P(2)^2*ydot(2)^2 - 2*P(2)*pm^2*x(2) + 2*P(2)*x(1)*ydot(1)*ydot(2) + 2*P(2)*x(2)*ydot(2)^2 + pm^2*x(1)^2 + pm^2*x(2)^2 - x(1)^2*ydot(1)^2 - 2*x(1)*x(2)*ydot(1)*ydot(2) - x(2)^2*ydot(2)^2)/(2*(P(1)*ydot(1)^3 + P(2)*ydot(2)^3 - x(1)*ydot(1)^3 - x(2)*ydot(2)^3 + P(1)*ydot(1)*ydot(2)^2 + P(2)*ydot(1)^2*ydot(2) + pm^2*x(1)*ydot(1) + pm^2*x(2)*ydot(2) - x(1)*ydot(1)*ydot(2)^2 - x(2)*ydot(1)^2*ydot(2) - P(1)*pm^2*ydot(1) - P(2)*pm^2*ydot(2)));
    a(1) = ((ydot(1)^2 + ydot(2)^2)^2 - ydot(1)^2 - ydot(2)^2);
    a(2) = (2*(ydot(1)^2 + ydot(2)^2)*(ydot(1)*(P(1) - x(1)) + ydot(2)*(P(2) - x(2))) - 2*ydot(2)*(P(2) - x(2)) - 2*ydot(1)*(P(1) - x(1)));
    a(3) = (ydot(1)*(P(1) - x(1)) + ydot(2)*(P(2) - x(2)))^2 - (P(1) - x(1))^2 - (P(2) - x(2))^2;
 
    if abs(a(1))<thresh
        if abs(a(2))<thresh
            stationaryPoints = [];
            orders = [];
        else
            stationaryPoints = roots(a(2:3));
            orders = 1;
        end
    else
         stationaryPointsInit = roots(a);
         [stationaryPoints,orders]  = sortStationaryPoints(stationaryPointsInit);
    end
        
    %check for degenerates, and group nearby stationary points if there are
    %two (I guess there usually should be two)
%     [stationaryPoints,orders]  = sortStationaryPoints(stationaryPointsInit);
    
end