function h0 = NSDpathICv3( SPorder, G, startPoint, ascFlag )
    %check that we've been given enough derivatives of phase:
    MinGsize = 2*SPorder-1;
    if length(G)<MinGsize
        error(fprintf('Current phase function has only %d terms, but PathFinder requires at least %d to compute ICs of ODE',length(G),MinGsize));
    end
    
    %h0 = cell(SPorder,1);
    for m=1:SPorder%
        h0{m}(1)=startPoint; %easy one

        if SPorder>=2            %almost as easy one
            AscDescDir=1;
            if nargin == 4
                if ascFlag
                    AscDescDir=-1;
                end
            end

            %these are all the possible paths satisfying the SD DE
            h0{m}(2)=(AscDescDir*1i*factorial(SPorder)/G{SPorder+1}(startPoint)).^(1/SPorder) * exp(1i*2*pi*(m)/SPorder);
        end

        for v=3:SPorder
%             uses a simplified version of Faa di Bruno's formula, ignoring
%             lower order derivatives of phase g at stationry points (as these are
%             zero).
               r=SPorder; %to be consistent with my notes' notation
               [partitionTypes, freq] = partitionCounter(v+r-2, r+1);
               numPartitionTypes = length(partitionTypes);
               totSum=0;
               for n=1:numPartitionTypes
                   summand=freq(n)*G{length(partitionTypes{n})+1}(startPoint);
                   for p=1:length(partitionTypes{n})
                       summand=summand*h0{m}(1+partitionTypes{n}(p));
                   end
                   totSum=totSum+summand;
               end
               denom = nchoosek(r+v-2,v-1) * G{r+1}(startPoint) * (h0{m}(2))^(r-1);
               h0{m}(v)=totSum/denom;
        end
    end

end