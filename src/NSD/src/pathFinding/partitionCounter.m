function [P, freq] = partitionCounter(N, minSize)
%Doesn't actually compute every partition of a set {1,...,N}, but instead
%returns the frequencies of each partition type. By 'partition type' I mean a
%partition consiting of some M subsets, the first of size m_1, the second of
%size m_2, the Mth of size m_M, such that the convention used is that the sequence is 
%monotonitcally increasing. The vector [m_1,...,m_M] represents the
%partition type.

%the variable minSize means that all partitions of size below this value
%are excluded, useful for my specific problem of FDB, as some of the lower
%derivatives will vanish at these stationary points.
if nargin==1
    minSize=1;
end
    
    %get all canonical forms:
    Q = tupSumCanon(N);
    
    %remove partitions with fewer sets, if requested
    if minSize>1
        pCount=1;
        if minSize>N
            error('minsize must be no greater than N');
        end
        for q=1:length(Q)
           if length(Q{q}) >= minSize
               P{pCount}=Q{q};
               pCount=pCount+1;
           end
        end
    else
        P=Q;
    end
    
    %now count the frequencies, using similar ideas to period orbits of
    %permutation groups from Bunce's algebra course
    totalTypes=length(P);
    freq=ones(1,totalTypes);
    for j=1:totalTypes
       choices=N;
       repeats=ones(1,N);
        for p=1:length(P{j})
            freq(j)=freq(j)*nchoosek(choices,P{j}(p));
            %keep track of repetitions in the sequence, as these falsely
            %suggest that there are more possible partitions than there are
            if p>1 & P{j}(p)==P{j}(p-1)
                repeats(P{j}(p)) = repeats(P{j}(p))+1;
            end
            choices=choices-P{j}(p);
        end
        for r=repeats
            freq(j)=freq(j)/factorial(r);
        end
    end
end

