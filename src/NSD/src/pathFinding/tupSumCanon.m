function tups = tupSumCanon(N)
    tupleIndex=1;
    tups={};
    
    %go
    tupMaker(1);
    
    %function within a function, done to exploit global nature of the
    %tuple, whilst using recursive nature of N nested loops
    function tupMaker(tupleEntryIndex)
        
        if ~isempty(tups)
            %there are two ways a tuple can be completed:
            if tupleFull()
                tups{tupleIndex+1}=tups{tupleIndex}(1:(tupleEntryIndex-2));
                tupleIndex=tupleIndex+1;
                return;
            elseif tupleEntryIndex>=N
                tups{tupleIndex}(tupleEntryIndex)=N-sum(tups{tupleIndex});
                tups{tupleIndex+1}=tups{tupleIndex}(1:(tupleEntryIndex-1));
                tupleIndex=tupleIndex+1;
                return;
            end
        end
        
        %otherwise do this:
        if tupleEntryIndex==1
            nEnd=N;
            nStart=1;
        else
            nEnd=(N-sum(tups{tupleIndex}));
            nStart=tups{tupleIndex}(tupleEntryIndex-1);
            %Nlim=tups{tupleIndex}(tupleEntryIndex-1);
        end
        for n=nStart:nEnd
            tups{tupleIndex}(tupleEntryIndex)=n;
            if tups{end}(1)==N
                %then we've reached the last pertubation
                return;
            end
            %now do the next entry:
            tupMaker(tupleEntryIndex+1);
        end
    end

    %test to see if tuple is 'full'
    function tf=tupleFull()
        if sum(tups{tupleIndex})==N
            tf=true;
        else
            tf=false;
        end
    end

end