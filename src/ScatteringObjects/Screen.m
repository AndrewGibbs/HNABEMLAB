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
        
        function drawCompliment(self,sMin,sMax)
            seg = [0 self.L];
            S = sort([sMin seg(1 == ((seg<sMax) .* (seg>sMin))) sMax]);
            
            v1 = self.vertices(1,:);
            v2 = self.vertices(2,:);
            d = (v2-v1);
            
            nodes = v1 + S.'*d;
            
            hold on;
            for n = 1:(length(S)/2)
                %L = norm(nodes(2*n,:)-nodes(2*n+1,:));
                s = linspace(0,self.L).';
                X = linspace(nodes(2*n-1,1),nodes(2*n,1));
                Y = linspace(nodes(2*n-1,2),nodes(2*n,2));
                plot(X,Y,'k','LineWidth',3);
                clear s X Y;
            end
            hold off;
            
        end
    end
end

