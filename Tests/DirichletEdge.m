classdef DirichletEdge < EdgeFunction
    
    properties
        
    end
    
    methods
        function self = DirichletEdge(uinc,component)
            %now deal with abstract function properties
            self.oscillator=@(s) uinc.oscillator(component.trace(s));
            self.supp=[0 component.L];
            self.suppWidth=component.L;
        end
        
        function [yOsc, y]=eval(self, s)
            [yOsc, y]=self.uinc.DirTrace(s,self.domain);
        end
    end
end