function thePath = NSDpathICv2( r, InOut, G, SP, ascFlag )
        
    AscDescDir=1; %
    if nargin == 5
        if ascFlag
            AscDescDir=-1;
        end
    end
    %some small number to move along SD path by:
    small=1E-5;
    
    %these are all the possible paths satisfying the SD DE
    allPaths=(AscDescDir*1i*factorial(r)/G{r+1}(SP)).^(1/r) * exp(1i*2*pi*(0:(r-1))/r);
    
    if r==1
        thePath=allPaths;
        return;
    end
    
    if AscDescDir==1
        if InOut==1
            somePaths=allPaths(real(allPaths)>=0);
        else
            somePaths=allPaths(real(allPaths)<=0);
        end
        %somePaths=allPaths;
    else
        if InOut==1
            somePaths=allPaths(sign(imag(SP))==-sign(imag(allPaths)));
%         else
%             somePaths=allPaths(real(allPaths)<0);
        end
    end
    
    
    %now choose path with max decay/growth
    
    thePath=somePaths(1);
    minVal=abs(exp(1i*G{1}((SP+small*thePath)^r)));
    for j=2:length(somePaths)
       newMinValMaybz=abs(exp(1i*(G{1}(SP)+(small*somePaths(j))^r)));
       if newMinValMaybz<minVal
            thePath=somePaths(j);
            minVal=newMinValMaybz;
       end
    end

end