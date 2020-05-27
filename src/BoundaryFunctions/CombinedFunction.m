classdef CombinedFunction < BoundaryFunction
    
    properties
        kwave
        uinc
        eta
    end
    
    methods 
        function self = CombinedFunction(uinc,domain,eta)
            
            edgeComponents = @(x) CombinedEdge(uinc,x,eta);
            
            self@BoundaryFunction(domain, edgeComponents);
            
            self.uinc=uinc;
            self.kwave=uinc.kwave;
            self.domain=domain;
            
        end

    end
    
end

