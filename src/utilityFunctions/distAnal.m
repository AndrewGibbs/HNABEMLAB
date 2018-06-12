function R = distAnal(s,t,deriv,sGEt,sEdge,tEdge)
    %alt inputs: (s,x,deriv,[],tEdge,[])
    if nargin == 5 %alternate set of inputs, maybe could be used for point source etc
        if isequal(size(t),[1,2])
            tEdge = sEdge;
            x = t;
            clear sEdge;
        else
            error('Please provide second argument as 2D vector');
        end
    else
        x = sEdge.trace(s);
    end
    t=t(:);
    if isequal(sEdge,tEdge) %same edge
        R = sEdge.distAnal(self,s,t,deriv,sGEt);
    else  % some condition about the sides not being too close together...
        y = tEdge.trace(t);
        xMy = [ (x(1)-y(:,1))  (x(2)-y(:,2)) ];
        dy =tEdge.dSv;
        R0 = sqrt( (x(1) - y(:,1) ).^2 + (x(2) - y(:,2) ).^2);
        switch deriv
            case 0 
                R = R0;
            case 1
                R = (dy*xMy.').'./R0;
            case 2
                R = (dy * ( [dy(1)./R0 dy(2)./R0] - [y(:,1)./(R0.^3)  y(:,2)./(R0.^3)] * sum(xMy(1,:)*dy(1) + xMy(2,:)*dy(2)) ).').';
            case 3
                error('Havent coded derivtives this high yet')
        end
        %alternative option if the sides meet at a vertex, or are
        %close together
    end
end