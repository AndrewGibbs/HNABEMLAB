classdef GeometricalOpticsFunction < BoundaryFunction
    %only valid for screens!
    
    properties
        uinc
        %phaseLinear
        sourceVsnormal
        dirConst
        illumSides = []
    end
    
    methods
        function self=GeometricalOpticsFunction(uinc,domain)
            self.uinc=uinc;
            self.domain=domain;
%             if ~isequal(uinc.phaseLinear,[])
%                 self.phaseLinear=[  uinc.phaseLinear*domain.dSv.'  uinc.phaseLinear*domain.P1.' ];
%             else
%                 self.phaseLinear=[];
%             end
            
            %now deal with abstract function properties
            self.oscillator=@(s) self.uinc.oscillator(self.domain.trace(s));
            self.supp=[zeros(domain.numSides,1) domain.L];
            self.suppWidth=self.domain.L;
            
%             self.a = self.supp(1);
%             self.b = self.supp(2);
            self.phaseMaxStationaryPointOrder = uinc.phaseMaxStationaryPointOrder;
            
            %+-1 to determine direction of reflection
            if isa(domain,'edge')
                self.dirConst = -uinc.sourceVsnormal(domain);
                self.illumSides = 1;
            else %on a polygon, zero on sides which can't 'see' inc wave
                for n=1:domain.numSides
                    pre_dirConst = uinc.sourceVsnormal(domain.side{n});
                    if pre_dirConst<0
                        self.dirConst(n) = 1;
                        self.illumSides = [self.illumSides n];
                    else
                        self.dirConst(n) = 0;
                    end
                    clear pre_dirConst;
                end
            end
        end
        
        function NeuOsc = eval(self, s, onSide)
            %get edge of onSide
            Edge = self.domain.getEdge(onSide);
            %compute result
            NeuOsc=2*self.dirConst(onSide)*self.uinc.NeuTrace(s,Edge);
        end
            
        function valNonOsc = evalNonOscAnal(self, s, onSide)
            
            Edge = self.domain.getEdge(onSide);
            
            [~,NeuNonOsc]=self.uinc.NeuTrace(s, Edge);
            %double it
            valNonOsc=2*self.dirConst(onSide)*NeuNonOsc;
        end
        
        function g = phaseAnal(self,s,deriv,onSide)
            s=s(:);
            
            Edge = self.domain.getEdge(onSide);
            
            switch deriv
                case 0 
                    g = self.uinc.phasePD(Edge.trace(s,onSide),0,0);
                case 1
                    g = Edge.dSv(1)*self.uinc.phasePD(self.domain.trace(s,onSide),1,0) ...
                         + Edge.dSv(2)*self.uinc.phasePD(Edge.trace(s,onSide),0,1);
                case 2
                    g =  Edge.dSv(1)^2*self.uinc.phasePD(Edge.trace(s,onSide),2,0) ...
                          + Edge.dSv(2)^2*self.uinc.phasePD(Edge.trace(s,onSide),0,2) ...
                          + 2*Edge.dSv(1)*Edge.dSv(2)*self.uinc.phasePD(Edge.trace(s),1,1);
                otherwise
                    error('Havent coded for such high partial derivatives yet');
            end
        end
       
       
    end
    
end