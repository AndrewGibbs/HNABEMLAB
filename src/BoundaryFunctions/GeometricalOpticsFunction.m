classdef GeometricalOpticsFunction < BoundaryFunction
    %only valid for screens!
    
    properties
        uinc
        phaseLinear
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
            
%             self.phase{1} = @(s) uinc.phasePD{1,1}(self.domain.trace(s));
%             self.phase{2} = @(s) domain.dSv(1)*uinc.phasePD{2,1}(self.domain.trace(s)) ...
%                                  + domain.dSv(2)*uinc.phasePD{1,2}(self.domain.trace(s));
%             if length(uinc.phasePD{1,1})>2
%                 self.phase{3} = @(s)  domain.dSv(1)^2*uinc.phasePD{3,1}(self.domain.trace(s)) ...
%                                     + domain.dSv(2)^2*uinc.phasePD{1,3}(self.domain.trace(s)) ...
%                                     + 2*domain.dSv(1)*domain.dSv(2)*uinc.phasePD{2,2}(self.domain.trace(s));
%             end
%             if length(uinc.phasePD{1,1})>3
%                 error('Havent coded for such high partial derivatives yet');
%             end
            self.a = self.supp(1);
            self.b = self.supp(2);
        end
        
        function [val, valNonOsc] = eval(self,s)
            % ** keeping second argument of this code so old stuff still
            % runs, should delete this once the new stuff has been
            % thoroughky tested
            %get Neumann data
            [Neu,NeuNonOsc]=self.uinc.NeuTrace(s,self.domain);
            %double it
            val=2*Neu;
            valNonOsc=2*NeuNonOsc;
        end
            
        function valNonOsc = evalNonOscAnal(self,s)
            [~,NeuNonOsc]=self.uinc.NeuTrace(s,self.domain);
            %double it
            valNonOsc=2*NeuNonOsc;
        end
        
        function g = phaseAnal(self,s,deriv,~)
            switch deriv
                case 0 
                    g = self.uinc.phasePD{1,1}(self.domain.trace(s));
                case 1
                    g = self.domain.dSv(1)*self.uinc.phasePD{2,1}(self.domain.trace(s)) ...
                         + self.domain.dSv(2)*self.uinc.phasePD{1,2}(self.domain.trace(s));
                case 2
                    g =  self.domain.dSv(1)^2*self.uinc.phasePD{3,1}(self.domain.trace(s)) ...
                          + self.domain.dSv(2)^2*self.uinc.phasePD{1,3}(self.domain.trace(s)) ...
                          + 2*self.domain.dSv(1)*self.domain.dSv(2)*self.uinc.phasePD{2,2}(self.domain.trace(s));
                otherwise
                    error('Havent coded for such high partial derivatives yet');
            end
        end
       
    end
    
end