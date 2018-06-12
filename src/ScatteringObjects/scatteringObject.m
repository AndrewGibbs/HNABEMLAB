classdef (Abstract) scatteringObject < handle
    %abstract class of scattering objects, such as screens, convex
    %polygons, etc
    
    properties
        L
        supp
        numSides
        vertices
    end
    
    methods
        val = trace(self,s);
        val = normal(self,s);
        val = dist(self,s,t);
        val = distAnal(self,s,t,deriv,sGEt);
        
        function Edge = getEdge(self,side)
            if isa(self,'edge')
                Edge = self;
            elseif isa(self,'polygon')
                Edge = self.side{side};
            end
        end
        
        function draw(self)
           plot(self.vertices(:,1),self.vertices(:,2)) ;
        end
        
        
    end
end

