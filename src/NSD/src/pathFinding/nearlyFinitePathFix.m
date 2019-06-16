function [h, dhdp, w] = nearlyFinitePathFix(nfParamPoint, nfCritPoint, criticalStartPoints, pathPower, IC, G, freq, N, COVsingularities, RelTol)

    %have taken special care to avoid including the endpoint twice
    
    Pend=nfParamPoint*2^(-1/pathPower); %truncate the path far from the critical point,
    %it doesn't matter where too much, as long as it's less than nfParamPoint

    %make the first bit of the path:
    [p, Wstart] = pathQuadFinite( Pend, COVsingularities, pathPower, freq, N );
    p_0pEnd = [0; p ;Pend];
    [~, Hstart] = ode45(@(t,y) NSDpathODE(t,y,pathPower-1,G, IC, false), p_0pEnd, IC, odeset('RelTol',RelTol) );   
    %p_0pEnd
    %now go from here to the nearby critical point, in a straight line:
    L=abs(Hstart(end,1)-nfCritPoint);
    dhEnddp=(nfCritPoint-Hstart(end,1))/L;
    [zEnd, Wend] = quad_gauss(N, 0, L);
    Hend = Hstart(end,1) + zEnd*dhEnddp;
    
    %piece this wonky path together:
    h = [Hstart(2:(end-1),1); Hend];
    w = [Wstart*exp(1i*freq*G{1}(IC(1))); Wend.*exp(1i*freq*G{1}(Hend))];% 
    if pathPower ==1 
        dhdp = 1i./G{2}(h); %back into the ODE 
    else
        dhdp = [Hstart(2:(end-1),2); ones(length(Wend),1)*dhEnddp];
    end
    
end