function [B,A] = getSymmetryIndices(Op,N,M)
%predicts repitions in collocation matrix and exploits this
    X = Op.getSymmetries();
    numComps = Op.domain.numComponents;
    DOFsPerComp = N/numComps;
    COLsPerComp = M/numComps;
    A = NaN(M,N);
    B = NaN(M,N);
    parentBlocks = [];
    parentInds = [];
    parentCount = 0;
    for n = 1:numComps
        for m = 1:numComps
            if ismember(X(m,n),parentBlocks)
               [aInds, bInds] = getFullsIndices(m,n);
               A(aInds,bInds) = repmat(parentInds{X(m,n)}.A,[DOFsPerComp 1]).';
               B(aInds,bInds) = repmat(parentInds{X(m,n)}.B.',[1 COLsPerComp]).';
            else
                parentBlocks = [parentBlocks X(m,n)];
                parentCount = parentCount + 1;
                [aInds, bInds] = getFullsIndices(m,n);
                parentInds{X(m,n)}.A = aInds;
                parentInds{X(m,n)}.B = bInds;
            end
        end
    end
    
    function [aInds, bInds] = getFullsIndices(m,n)
        bInds = DOFsPerComp*(n-1)+1 : DOFsPerComp*n;
        aInds = COLsPerComp*(m-1)+1 : COLsPerComp*m;
    end
end