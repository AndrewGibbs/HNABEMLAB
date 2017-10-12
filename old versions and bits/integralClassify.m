function type=integralClassify(X,Y, nearParam)
%assumes integrals have already been split into subintegrals
    if nargin==2
        nearParam=0.15;
    end
    if isequal(X,Y)
        type='diagSing';
    elseif X(1)==Y(2) || Y(1)==X(2)
        type='cornerSing';
    elseif min(abs(X(1)-Y(2)), abs(Y(1)-X(2)))<nearParam;
        type='nearSing'; 
    else
        type='smooth';
    end
end