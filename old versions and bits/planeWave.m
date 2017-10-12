classdef planeWave < incField
    %plane wave incident field
    
    properties
        d
    end
    
    methods
        
        function self=planeWave(kwave,direction)
           self.kwave=kwave;
           self.d=direction;
        end
        
        function val=DirTrace(self,s,boundary)
            %standard def of plane wave
            val=exp(1i*self.kwave*boundary.trace(s)*self.d.');
        end
        
        function val=NeuTrace(self,s,boundary)
            %derivative of plane wave, product with normal vector
            val=1i*self.kwave*boundary.normal(s)*self.d.'.*self.DirTrace(s,boundary);
        end
    end
    
end

