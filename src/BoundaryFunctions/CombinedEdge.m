classdef CombinedEdge < EdgeFunction
    
    properties
        uinc
        eta
    end
    
    methods
        function self = CombinedEdge(uinc,component,eta)
            self.uinc = uinc;
            %now deal with abstract function properties
            self.oscillator=@(s) uinc.oscillator(component.trace(s));
            self.supp=[0 component.L];
            self.suppWidth=component.L;
            self.domain = component;
            self.eta = eta;
        end
        
        function [yOsc, y] = eval(self, s)
            [yOsc_d, y_d] = self.uinc.DirTrace(s,self.domain);
            [yOsc_n, y_n] = self.uinc.NeuTrace(s,self.domain);
            yOsc = yOsc_n - 1i*self.eta*yOsc_d;
            y    = y_n    - 1i*self.eta*y_d;
        end
    end
end