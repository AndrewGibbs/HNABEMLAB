function [gStationaryPoints, gSingularities, gSPorders,gPoleOrders] = getStationaryPoints(a,b,rectRad,...
                                                                     analytic, G, RectTol, Npts , visuals,...
                                                                     settleRad, intThresh)
                
    %Requires at least up to G{3}:=g"(x)
     if length(G)<3
            error(['Require up to second derivative of phase non-analytic functions' ...
                'OR an anlytic phase function, to automatically detect stationary points']);
     end
     
     %This one should work better, although the '2' is quite arbitrary.
    if isempty(rectRad)
        if isinf(a+b)
            rectRad=settleRad;
        else
            rectRad=.75*(b-a);
        end
    end
    
    %now check that phase isn't a constant function, as this could mess
    %things up if g'(z)=0 on compact subset of rectangle
    if constantPhase(a,b,G,Npts,intThresh)
        gStationaryPoints = [];
        gSingularities = [];
        gSPorders = []; 
        gPoleOrders = [];
        return; %early exit
    end
    
    %scale quad points by length of rectangle:
    Nrect=Npts*max(1,ceil(rectRad));
    
    if isinf(a+b)
       a=0; b=0; %centre rectangle at origin, as there's no other point of reference
    end
    
    initRect=[a-rectRad-rectRad*1i  b+rectRad-rectRad*1i  b+rectRad+rectRad*1i  a-rectRad+rectRad*1i];  

    if analytic %can do a quicker verison of the zeros search:
        [gStationaryPoints, gSPorders] = findZerosRect( G{2}, G{3}, initRect, RectTol, intThresh, Nrect, visuals );
        gSingularities=[];
        gPoleOrders=[];
    else %need ot zoom in further into the rectangle:
        %now find all stationary points inside of this rectangle
        [gStationaryPoints, gSingularities, gSPorders, gPoleOrders] = findZerosSingsRect( G{2}, G{3}, initRect, RectTol, Nrect , visuals);        
    end

end