classdef DirichletFunction < BoundaryFunction
    %data for standard combined formulation, (typically denoted f_{k,\eta}
    
    properties
        kwave
        uinc
        phaseLinear
    end
    
    methods 
        function self=DirichletFunction(uinc,domain)
            self.uinc=uinc;
            self.kwave=uinc.kwave;
            self.domain=domain;
            
            %now deal with abstract function properties
            self.oscillator=@(s) self.uinc.oscillator(self.domain.trace(s));
            self.supp=[zeros(domain.numSides,1) domain.L];
            self.suppWidth=self.domain.L;

%             if ~isequal(uinc.phaseLinear,[])
%                 self.phaseLinear=[self.domain.dSv*uinc.phaseLinear.'  self.domain.P1*uinc.phaseLinear.'];
%             end
        end
        
        function [yOsc, y]=eval(self,s, sSide)
            if nargin ==2 || isa(self.domain,'edge')
                [yOsc, y]=self.uinc.DirTrace(s,self.domain);
            else
                [yOsc, y]=self.uinc.DirTrace(s,self.domain.side{sSide});
            end
        end

    end
    
end

