classdef DirichletEdge < EdgeFunction
    
    properties
        uinc
    end
    
    methods
        function self = DirichletEdge(uinc,component)
            self.uinc = uinc;
            %now deal with abstract function properties
            self.oscillator=@(s) uinc.oscillator(component.trace(s));
            self.supp=[0 component.L];
            self.suppWidth=component.L;
            self.domain = component;
        end
        
        function [yOsc, y]=eval(self, s)
            [yOsc, y]=self.uinc.DirTrace(s,self.domain);
        end
    end
end