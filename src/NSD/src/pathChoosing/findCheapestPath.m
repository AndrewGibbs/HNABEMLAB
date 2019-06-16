function [path, cost] = findCheapestPath( P, g,  freq, Npts, m, pathPowers, FPfullIndices)
%looks at start and 'end' points of each SD path computed, and finds the
%shortest path from a to b

    Plength=size(P); Plength=Plength(1);
    m=Plength;
    A=inf(Plength);
    small=1E-14; %this amount needs to be small, but not too small.
    %All edges need to have some length, but if that length is negligable,
    % some silly connected paths can occur, especially when there are finite paths,
    % as an SD path can end at the beginning at another
    
    for j=1:m
       for ell=(j+1):m
           if P(j,1)==P(ell,1) %start points are the same
              A(ell,j)=small; %need this bodge - matlab only connects non-zero entries
              
              %next elseif contains some confusing stuff. The point is
              %to map the path indices to and from the critical point
              %indices, the latter of which is fewer than the former
           elseif ismember(ell,FPfullIndices{j}) || ismember(j,FPfullIndices{ell})
                  A(ell,j)=small;
           else   %endpoints are possibly 'connected'
              A(ell,j)=pathCost(P(j,2),P(ell,2), g,  freq, Npts);  
           end             
           A(j,ell)=A(ell,j); %symmetric obvz
       end
    end
    
    %now compute shortest (cheapest) path from a to b
    G = graph(A.','upper');
    
    cost=inf;
    for startPath=1:pathPowers(1)
        for endPath=(m-pathPowers(end)+1):m
            [path_, cost_] =shortestpath(G,1,m);
            if cost_<cost
                path=path_;
                cost=cost_;
            end
        end
    end
    
end