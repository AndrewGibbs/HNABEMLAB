function h0Derivs = NSDpathICvAllTerms( SPorder, G, startPoint, ascFlag )

    NumTerms = length(G)-SPorder;
    
    for m=1:SPorder%
        h0Derivs{m}(1)=startPoint; %easy one

        if NumTerms>=2            %almost as easy one
            AscDescDir=1;
            if nargin == 4
                if ascFlag
                    AscDescDir=-1;
                end
            end

            %these are all the possible paths satisfying the SD DE
            h0Derivs{m}(2)=(AscDescDir*1i*factorial(SPorder)/G{SPorder+1}(startPoint)).^(1/SPorder) * exp(1i*2*pi*(m)/SPorder);
        end

        for v=3:NumTerms
               deriv = v-1;
               r=SPorder; %g^(r)(\xi) \neq 0
               [partitionTypes, freq] = partitionCounter(deriv+r-1, r);
               numPartitionTypes = length(partitionTypes);
               totSum=0;
               denomCount = 0;
               for n=1:numPartitionTypes
                   if ismember(deriv,partitionTypes{n}) %the term we want to divide through by
                       denom = freq(n) * G{r+1}(startPoint) * (h0Derivs{m}(2))^(r-1);
                       denomCount = denomCount + 1;
                       if denomCount>1
                           error('too many terms to divide by');
                       end
                   else %all the other terms
                       summand = freq(n)*G{length(partitionTypes{n})+1}(startPoint);
                       for p=1:length(partitionTypes{n})
                           summand = summand * h0Derivs{m}(1+partitionTypes{n}(p));
                       end
                       totSum = totSum + summand;
                   end
               end
               %denom = nchoosek(r+v-2,v-1) * G{r+1}(startPoint) * (h0Derivs{m}(2))^(r-1);
               h0Derivs{m}(v) = -totSum/denom;
        end
    end

end