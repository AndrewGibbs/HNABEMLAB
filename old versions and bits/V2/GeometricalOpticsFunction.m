classdef GeometricalOpticsFunction < BoundaryFunction
    %only valid for screens!
    
    properties
        uinc
    end
    
    methods
        function self=GeometricalOpticsFunction(uinc,domain)
            self.uinc=uinc;
            self.domain=domain;
            
            %now deal with abstract function properties
            self.oscillator=@(s) self.uinc.oscillator(self.domain.trace(s));
            self.supp=[0 domain.L];
            self.suppWidth=self.domain.L;
        end
        
        function val=eval(self,s)
            val=2*self.uinc.NeuTrace(s,self.domain);
        end
    end
    
end