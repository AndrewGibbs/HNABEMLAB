classdef DirichletFunction < BoundaryFunction
    %data for standard combined formulation, (typically denoted f_{k,\eta}
    
    properties
        kwave
        uinc
    end
    
    methods 
        function self=DirichletFunction(uinc,domain)
            self.uinc=uinc;
            self.kwave=uinc.kwave;
            self.domain=domain;
            
            %now deal with abstract function properties
            self.oscillator=@(s) self.uinc.oscillator(self.domain.trace(s));
            self.supp=[0 domain.L];
            self.suppWidth=self.domain.L;
        end
        
        function [yOsc]=eval(self,s)
            yOsc=self.uinc.DirTrace(s,self.domain);
        end
    end
    
end

