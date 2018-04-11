function GOA=GeometricalOpticsApprox(uinc,domain)
    %decides on appropriate class for GOA given uinc
    if isa(uinc,'BoundaryIntegral')
        GOA=2*uinc.NeuTrace(domain);
    elseif isa(uinc,'waveR2')
        GOA=GeometricalOpticsFunction(uinc,domain);
    else
        error('Inc field must be waveR2 or BoundaryIntegral type');
    end
end