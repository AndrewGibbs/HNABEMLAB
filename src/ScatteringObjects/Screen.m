classdef Screen < PolygonalScatteringObject
    %contains all information and methods needed about a polygon
    
    properties
        L %some measure of total length of boundary
        nv
    end
    
    methods
        function self = Screen(vertices)
            
            self.component(1) = edge(vertices);
            
            self.numComponents = 1;
            
            
            self.vertices = vertices;
            
            self.L = self.component(1).L;
            self.nv = self.component(1).nv;
            
            %now make matrix of internal angles, easy in this case:
            self.internalAngle = pi;
            
            self.Lipschitz = false;
            
        end
        
    end
end

