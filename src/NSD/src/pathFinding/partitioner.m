function B = partitioner(N, minSize)
%returns all partitions of the set \PI={1,...,N}
    if nargin==1
        minSize=1;
    end
    PI=1:N;
    indexTuples=tupSum(N);
    tupCount=length(indexTuples);
    
    partCount=0;
    for t=1:tupCount
        if length(indexTuples{t})>=minSize
            partCount=partCount+1;
            indices=0;
            Tcount=0;
            for T=indexTuples{t}
                indices=(indices(end)+1):indices(end)+T;
                Tcount=Tcount+1;
                B{partCount}{Tcount}=PI(indices);
            end
        end
    end
    
end