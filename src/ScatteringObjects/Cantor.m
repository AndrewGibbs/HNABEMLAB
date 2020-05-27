function s = Cantor(J,L)
%gives set of endpoints of prefractal cantor set on [0,L]
    if nargin == 1
        L = 1;
    end
    
    s = [0 1];
    removeBit = [1/3 2/3];
    for j=1:J
        for n=1:(length(s)/2)
            %get length of segment
            Ln = s(2*n) - s(2*n-1);
            sNew{n} = s(2*n-1) + removeBit*Ln;
        end
        for n = 1:length(sNew)
            s = sort([s sNew{n}]);
        end
    end
    
    s = L*s;
end

