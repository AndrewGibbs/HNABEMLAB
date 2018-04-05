classdef GeometricalOpticsFunction < BoundaryFunction
    %only valid for screens!
    
    properties
        uinc
        phaseLinear
        phase
    end
    
    methods
        function self=GeometricalOpticsFunction(uinc,domain)
            self.uinc=uinc;
            self.domain=domain;
            if ~isequal(uinc.phaseLinear,[])
                self.phaseLinear=[  uinc.phaseLinear*domain.dSv.'  uinc.phaseLinear*domain.P1.' ];
            else
                self.phaseLinear=[];
            end
            
            %now deal with abstract function properties
            self.oscillator=@(s) self.uinc.oscillator(self.domain.trace(s));
            self.supp=[0 domain.L];
            self.suppWidth=self.domain.L;
            
            phase = G{@(s) uinc.phase(trace(s)), };
        end
        
        function [val, valNonOsc]=eval(self,s)
            %get Neumann data
            [Neu, NeuNonOsc]=self.uinc.NeuTrace(s,self.domain);
            %double it
            val=2*Neu; 
            valNonOsc=2*NeuNonOsc;
        end
       
    end
    
end