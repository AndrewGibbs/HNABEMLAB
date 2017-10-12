classdef beamSource < incField
    %plane wave incident field
    
    properties
        sourceObj
        density
    end
    
    methods
        
        function self=planeWave(kwave,sourceObj_,density_)
           self.kwave=kwave_;
           self.sourceObj=sourceObj_;
           self.density=density_;
           self.asOperator=SingleLayer(self.kwave);
        end
        
        function val=DirTrace(self,s,boundary)
            %standard def of plane wave
            
        end
        
        function val=NeuTrace(self,s,boundary)
            %derivative of plane wave, product with normal vector
            val=1i*self.kwave*boundary.normal(s)*self.d.'.*self.DirTrace(s,boundary);
        end
    end
    
end

