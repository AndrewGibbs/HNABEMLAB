classdef StandardCombinedData < BoundaryFunction
    %data for standard combined formulation, (typically denoted f_{k,\eta}
    
    properties
        kwave
        eta     %coupling parameter
        uinc
    end
    
    methods 
        function self=StandardCombinedData(uinc,domain,eta)
            if nargin<=2
                warning('No coupling parameter specified, set to wavenumber.');
                self.eta=uinc.kwave;
            else
                self.eta=eta;
            end
            self.uinc=uinc;
            self.kwave=uinc.kwave;
            self.domain=domain;
            
            %now deal with abstract function properties
            self.oscillator=@(s) self.uinc.oscillator(self.domain.trace(s));
            self.supp=[0 domain.L];
            self.suppWidth=self.domain.L;
        end
        
        function [yOsc]=eval(self,s)
            yOsc=self.uinc.NeuTrace(s,self.domain)-1i*self.eta*self.uinc.DirTrace(s,self.domain);
        end
    end
    
end

