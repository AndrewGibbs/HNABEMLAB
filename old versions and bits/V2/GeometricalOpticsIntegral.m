classdef GeometricalOpticsIntegral < BoundaryIntegral
    %only valid for screens!
    
    properties
        uinc
    end
    
    methods
        function self=GeometricalOpticsIntegral(uinc,domain)
            %here uinc is also a boundary integral on R2
            self@BoundaryIntegral(kernel,boundaryFn,boundary);
            self.uinc=uinc;
            self.domain=uinc.domain;
            self.domain{length(self.domain)+1}=domain; %trace taken here
            self.kernel=2*uinc.kernel; %WILL CAUSE PROBLEMS
            self.Fn=uinc.Fn;
            
            %now deal with abstract function properties
            self.oscillator=@(s) self.uinc.oscillator(self.domain.trace(s));
            self.supp=[0 domain.L];
            self.suppWidth=self.domain.L;
        end
        
%         function val=eval(self,s)
%             val=2*self.uinc.NeuTrace(s,self.domain);
%         end
    end
    
end