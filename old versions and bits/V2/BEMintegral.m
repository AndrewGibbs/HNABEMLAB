classdef (Abstract) BEMintegral
    %generic integrator compnent which is intentded to be swapped in and
    %out
    
    properties
        domain
        kernel
        oscillator
        type
        toR2    %flag which says if this maps to function in R2 or boundary
    end
    
    methods
    end
    
end

