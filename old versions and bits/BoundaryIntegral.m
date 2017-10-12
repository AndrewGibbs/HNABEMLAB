classdef BoundaryIntegral
    %general class which 
    
    properties
        domain
        kernel
        boundaryFn
    end
    
    methods 
        function self=BoundaryIntegral()
            
        end
        
        function val=eval(self,pt,integrator)
            switch length(domain)
                case 1
                    
                case 2
                
            end
        end
    end
    
end

