classdef GeometricalOpticsEdge < EdgeFunction
    %only valid for screens!
    
    properties
        uinc
        sourceVsnormal
        dirConst = []
        nodes
        weights
    end
    
    methods
        function self=GeometricalOpticsEdge(uinc,edge)
            self.uinc=uinc;
            self.domain=edge;
            
            %now deal with abstract function properties
            self.oscillator=@(s) self.uinc.oscillator(self.domain.trace(s));
            self.supp=[zeros(edge.numSides,1) edge.L];
            self.suppWidth=self.domain.L;
            self.a = self.supp(1);
            self.b = self.supp(2);
            
            self.phaseMaxStationaryPointOrder = uinc.phaseMaxStationaryPointOrder;
            self.meshEl = wholeMeshSide(edge.L);
            
            [self.nodes, self.weights] = gauss_quad_wave_split2(self.a, self.b, self.nodesPerWavelength,  uinc.kwave, edge.L );
        end
        
        function val = eval(self, s)
            if nargin == 2
                onSide = 1;
            end
              
            %get edge of onSide
            Edge = self.domain;
            %compute result
            val = 2*self.dirConst*self.uinc.NeuTrace(s,Edge);
        end
            
        function valNonOsc = evalNonOscAnal(self, s, ~)
            
            Edge = self.domain;
            
            [~,NeuNonOsc]=self.uinc.NeuTrace(s, Edge);
            %double it
            valNonOsc=2*self.dirConst*NeuNonOsc;
        end
        
        function g = phaseAnal(self,s,deriv,onSide)
            s=s(:);
            
            Edge = self.domain;
            
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