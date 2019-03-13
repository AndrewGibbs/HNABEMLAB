classdef Screen < PolygonalScatteringObject
    %contains all information and methods needed about a polygon
    
    properties
        L %some measure of total length of boundary
        nv
    end
    
    methods
        function self = Screen(vertices)
            [sizeY, sizeX]= size(vertices);
            if sizeX ~= 2 || sizeY ~= 2
                error('First input must be 2x2 vector');
            end
            
            self.numComponents = 1;
            
            self.component{1} = edge(vertices);
            
            self.vertices = vertices;
            
            self.L = self.component{1}.L;
            self.nv = self.component{1}.nv;
            
            %now make matrix of internal angles, easy in this case:
            self.internalAngle = pi;
            
            self.Lipschitz = false;
            
        end
        
    end
end

