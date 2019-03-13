classdef DirichletFunction < BoundaryFunction
    
    properties
        kwave
        uinc
    end
    
    methods 
        function self=DirichletFunction(uinc,domain)
            
            edgeComponents = @(x) DirichletEdge(uinc,x);
            
            self@BoundaryFunction(domain, edgeComponents);
            
            self.uinc=uinc;
            self.kwave=uinc.kwave;
            self.domain=domain;
            
        end

    end
    
end

