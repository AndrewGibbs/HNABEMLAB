function orders = stationaryPointsGetOrders(sp,G)
%determines the order of a stationary point
thresh = 1e-7;
orders = zeros(size(sp));
    for m = 1:length(sp)
        x = sp(m);
        for n = 2:length(G)
            if abs(G{n}(x))<thresh
                orders(m) = orders(m) + 1;
            else
                break;
            end
        end
    end

end

