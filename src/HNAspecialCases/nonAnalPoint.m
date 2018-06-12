function z = nonAnalPoint(x,Edge)
%computes the points at which the analytic extension of the distance
%function is not analytic
    P=Edge.P1; P1=P(1); P2=P(2);
    ydot=Edge.dSv; yd1=ydot(1); yd2=ydot(2);
    x1 = x(1); x2=x(2);
     
    z(1) = -(P1 - P2*1i - x1 + x2*1i)/(yd1 - yd2*1i);
    z(2) =  -(P1*1i - P2 - x1*1i + x2)/(yd1*1i - yd2);
end

