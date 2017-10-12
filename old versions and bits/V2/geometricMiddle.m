function b = geometricMiddle( Q, J, Q_width, J_width, near)
        lamda=( max([Q J])-min([Q J]))/( Q_width + J_width );
        %nearly touching
        if Q(2)<J(1)  && abs(J(1)-Q(2))<=near
            b=Q(1)+Q_width*lamda;
        elseif abs(Q(1)-J(2))<=near
            b=J(1)+J_width*lamda;
        else
            error('Cannot find geometric middle, or not near');
        end

end

