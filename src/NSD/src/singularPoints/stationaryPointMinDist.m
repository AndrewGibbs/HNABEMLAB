function p= stationaryPointMinDist(p,r,g,h01, h02, thresh)
%Netwon iteration along SD path to find point closest to a stationary point
%of ANOTHER SD path

    while err < thresh
        p = p - (abs(g(h01)-g(h02))^2 + p^(2*r))/(2*r*p^(2*r-1));
    end
end

