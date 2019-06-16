function branchPoints = getIntegrandBranchPoints(gSPorders,fSPorders, intThresh)
%determine non-singular branch points
    branchPoints=[];
    allSPorders=[gSPorders fSPorders];
    for j=1:length(allSPorders)
        if ~isNearlyInt( allSPorders(j), intThresh )
            branchPoints=[branchPoints allSPorders(j)];
        end
    end
end

