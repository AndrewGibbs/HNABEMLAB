function out=DirichletData(uinc,domain)
    %decides on appropriate class for GOA given uinc
    if isa(uinc,'BoundaryIntegral')
        out=uinc.DirTrace(domain);
    elseif isa(uinc,'waveR2')
        out=DirichletFunction(uinc,domain);
    else
        error('Inc field must be waveR2 or BoundaryIntegral type');
    end
end