function [stationaryPoints, orders, badPoints] = symbolicStationaryPointsCorner(colDistCorner, suppDistCorner, fun, internalAngle, funSide, phaseCorner)
%looks at types of phase and determines if stationary points can be
%determined symbolically

    if isa(fun,'GeometricalOpticsFunction')
        if isa(fun.uinc,'planeWave')
%             C = cos(internalAngle);
%             a = colDistCorner;
%             b = suppDistCorner;
%             d = fun.domain.side{funSide}.dSv*(fun.uinc.d.');
            
            %SP_flip1 = roots([1-d^2  -2*a*C*(1-d^2)  a^2*(d^2-C^2)+b*(1-d^2)-2*a*C*b*(1-d^2)]);
            d1 = fun.uinc.d(1); d2 = fun.uinc.d(2);
            cosTheta = cos(internalAngle) ;
            yd1 = fun.domain.side{funSide}.dSv(1);
            yd2 = fun.domain.side{funSide}.dSv(2);
            sx_dist = colDistCorner;
            ty_dist = suppDistCorner;
            if (d1*yd1 + d2*yd2)==1
                stationaryPoints = [];
                orders = [];
                return;
            else
            %elseif close2a
                SP_flip1(1) = (d1*sx_dist*yd1*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) - cosTheta*sx_dist + d2*sx_dist*yd2*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) + cosTheta*d1^2*sx_dist*yd1^2 + cosTheta*d2^2*sx_dist*yd2^2 + 2*cosTheta*d1*d2*sx_dist*yd1*yd2)/(d1^2*yd1^2 + 2*d1*d2*yd1*yd2 + d2^2*yd2^2 - 1) - ty_dist;
                SP_flip1(2) = - ty_dist - (cosTheta*sx_dist + d1*sx_dist*yd1*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) + d2*sx_dist*yd2*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) - cosTheta*d1^2*sx_dist*yd1^2 - cosTheta*d2^2*sx_dist*yd2^2 - 2*cosTheta*d1*d2*sx_dist*yd1*yd2)/(d1^2*yd1^2 + 2*d1*d2*yd1*yd2 + d2^2*yd2^2 - 1);
%             else
%                 SP_flip1(1) =  (d1*sx_dist*yd1*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) - cosTheta*sx_dist + d2*sx_dist*yd2*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) + cosTheta*d1^2*sx_dist*yd1^2 + cosTheta*d2^2*sx_dist*yd2^2 + 2*cosTheta*d1*d2*sx_dist*yd1*yd2)/(d1^2*yd1^2 + 2*d1*d2*yd1*yd2 + d2^2*yd2^2 - 1) - ty_dist;
%                 SP_flip1(2) =  - ty_dist - (cosTheta*sx_dist + d1*sx_dist*yd1*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) + d2*sx_dist*yd2*((cosTheta - 1)*(cosTheta + 1)*(d1*yd1 + d2*yd2 - 1)*(d1*yd1 + d2*yd2 + 1))^(1/2) - cosTheta*d1^2*sx_dist*yd1^2 - cosTheta*d2^2*sx_dist*yd2^2 - 2*cosTheta*d1*d2*sx_dist*yd1*yd2)/(d1^2*yd1^2 + 2*d1*d2*yd1*yd2 + d2^2*yd2^2 - 1);
            end
            %now sort these to make sure they're disjoint
            [SP_flip2, ~] = sortStationaryPoints(SP_flip1,1E-6);
             orders1 = stationaryPointsGetOrders(SP_flip2, phaseCorner);
            %delete the dodgy wrong ones:
            [stationaryPoints, orders] = deleteFalseStationaryPoints(phaseCorner{2}, SP_flip2, orders1);
                
        elseif isa(fun.uinc,'pointSource')
            stationaryPoints = NaN;
            orders = NaN;
        end
    elseif isa(fun,'baseFnHNA')
        stationaryPoints = [];
        orders = [];
    end
    
    badPoints = roots([1  2*suppDistCorner-2*colDistCorner*cos(internalAngle) suppDistCorner^2+colDistCorner^2-2*suppDistCorner*colDistCorner*cos(internalAngle)]);
    
%     %now sort these to make sure they're disjoint
%     [SP_flip2, orders2] = sortStationaryPoints(SP_flip1,1E-6);
%     %delete the dodgy wrong ones:
%     [SP_flip3, orders3] = deleteFalseStationaryPoints(phaseCorner{2}, SP_flip2, orders2);
end