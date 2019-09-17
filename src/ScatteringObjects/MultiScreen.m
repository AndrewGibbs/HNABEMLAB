classdef MultiScreen < PolygonalScatteringObject
    %contains all information and methods needed about a polygon
    
    properties
        L %some measure of total length of boundary
        nv
        segSplits
    end
    
    methods
        function self = MultiScreen(vertices, segSplits)
            
            [sizeY, sizeX]= size(vertices);
            if sizeX ~= 2 || sizeY ~= 2
                error('First input must be 2x2 vector');
            end
            
            self.numComponents = length(segSplits)/2;
            if mod(length(self.numComponents),1) ~ 0
               error('length of segSplits must be even');
            end
            
            self.L=zeros(self.numComponents,1);
            
            %create each component of screen as an 'edge'
            v1 = vertices(1,:);
            v2 = vertices(2,:);
            d = (v2-v1);
            segSplitsSorted = sort(segSplits-min(segSplits))/max(segSplits);
            
            self.segSplits = segSplitsSorted;
            
            self.L = 0;
            for j = 1:self.numComponents
                v1j = v1 + d*segSplitsSorted(2*j-1);
                v2j = v1 + d*segSplitsSorted(2*j);
                self.component(j) = edge([v1j; v2j]);
                self.L = self.L + self.component(j).L;
                %self.component{j}.nv = - self.component{j}.nv;
            end
            
            %now make matrix of internal angles, easy in this case:
            self.internalAngle = pi*ones(self.numComponents);
            
            self.Lipschitz = false;
            self.nv = self.component(1).nv;
            self.vertices = vertices;
        end
        
        function val = trace(self,s,ofSide)
            val = self.component(ofSide).trace(s);
        end
        
        function drawCompliment(self,sMin,sMax)
            S = sort([sMin self.segSplits(1 == ((self.segSplits<sMax) .* (self.segSplits>sMin))) sMax]);
            
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

