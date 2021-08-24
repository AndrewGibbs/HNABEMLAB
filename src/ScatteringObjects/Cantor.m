function s = Cantor(J,r)
%gives set of endpoints of prefractal cantor set on [0,L]   
L=1;
    if nargin ==1
        r = 1/3; 
    end
%     x = (1-r)/2;
    removeBit = L*[r 1-r];
    
    s = [0 1];
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

