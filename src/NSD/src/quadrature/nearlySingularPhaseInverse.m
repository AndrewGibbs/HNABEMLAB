function hNearSing = nearlySingularPhaseInverse(CP, SPs, threshold)

    hNearSing = [];
    
    for n = 1:length(SPs)
        if abs(SPs(n)-CP) < threshold && SPs(n)~=CP
            hNearSing = [hNearSing singularity(SPs(n), '??', @(r) abs(SPs(n)-CP)) ];
        end
    end
end