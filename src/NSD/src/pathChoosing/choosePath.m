function [X, W] = choosePath(a,b,P, G, freq, N, numPathsSD, pathPowers, visuals, X_, W_, FPfullIndices, R, ainf, binf)
       
    if ainf
       A=a*R;
       P=[A A; P];
       pathPowers=[1 pathPowers];
       FPfullIndices_{1}=[];
       for j=1:length(FPfullIndices)
           %shift LHS and RHS along by one
            FPfullIndices_{j+1}=FPfullIndices{j}+1;
       end
       FPfullIndices=FPfullIndices_;
    end
    if binf
       B=b*R;
       P=[P; B B];
       pathPowers=[pathPowers 1];
       FPfullIndices{length(FPfullIndices)+1}=[];
%        for j=1:length(FPfullIndices)
%            %shift RHS along by one
%             FPfullIndices{j}=FPfullIndices{j}+1;
%        end
    end
    
    [pathOrder, pathCost]=findCheapestPath( P, G{1}, freq, N, numPathsSD, pathPowers, FPfullIndices );%findFullPath( P, G{1}, freq, 1E-10, N, numPathsSD );
        
        if pathCost>.1
            warning(sprintf('Cost of truncating SD path is %e, may not be optimal path, %d %d %d %d',pathCost, pathOrder(1),pathOrder(2),pathOrder(3),pathOrder(4)));
        end

    %exlude the endpoints from the path order, as they were only artificial
    if ainf
       pathOrder=pathOrder(2:end)-1;
    end
    if binf
       pathOrder=pathOrder(1:(end-1));
    end
    
    X=[];   W=[]; xv=[]; yv=[];
    %if path starts at infinity, first SD contour should be negatively
    %weighted
    if ainf
        inOut=-1;
    else
        inOut=1;
    end
    
    for fullIndex=pathOrder
        W=[W; inOut*W_{fullIndex}];
        X=[X; X_{fullIndex};];       
        if inOut==1
            xv=[ xv; real(X_{fullIndex})];
            yv=[yv; imag(X_{fullIndex})];
        else
            xv=[ xv; real(flipud(X_{fullIndex}))];
            yv=[ yv; imag(flipud(X_{fullIndex}))];            
        end
        inOut=inOut*-1;
    end
    if visuals
        for fullIndex=1:numPathsSD
            if ismember(fullIndex,pathOrder)
                plot(X_{fullIndex},'xm'); 
            else
                plot(X_{fullIndex},'oc'); 
            end
        end
        %plot([a b],[0 0], 'b', 'LineWidth',2);
        if ~isempty(R)
            t=linspace(0,2*pi,1000);
            plot(R*exp(1i*t),'r--');
            plot(P(:,1),'ko')
        end
    end
    
end

