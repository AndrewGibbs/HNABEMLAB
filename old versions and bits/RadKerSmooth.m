classdef RadKerSmooth < BEMintegral
    properties
    end
    
    methods
        function RadKerSmooth(domain, kernel)
            %for domain and disjoint integrals
            type='smooth';
            self.domain=domain;
            self.kernel=kernel; %must be radial kernel
        end
    end
    
end

