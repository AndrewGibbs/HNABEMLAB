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
        
%         function draw(self)
%             for s = segSplits
%                 self.vertices(:,1)
%                 hold on;
%             end
%         end
        
%         function R = distAnal(self,s,t,deriv,sGEt,sSide,tSide)
%             s=s(:); t=t(:);
%             if sSide == tSide
%                 R = self.component(sSide).distAnal(s,t,deriv,sGEt);
%             else  % some condition about the sides not being too close together...
%                 x = self.component(sSide).trace(s);
%                 y = self.component(tSide).trace(t);
%                 yMx = -[ (x(1)-y(:,1))  (x(2)-y(:,2)) ];
%                 dy = self.component(tSide).dSv;
%                 R0 = sqrt( (x(1) - y(:,1) ).^2 + (x(2) - y(:,2) ).^2);
%                 switch deriv
%                     case 0 
%                         R = R0;
%                     case 1
%                         R = (dy*(yMx.')).'./R0;
%                     case 2
%                         %account for 2nd derivative bug:
%                         if length(t) == 1
%                             yMx1 = yMx(1);
%                             yMx2 = yMx(2);
%                         else
%                             yMx1 = yMx(1,:);
%                             yMx2 = yMx(2,:);
%                         end
%                         R = - dy(1)*((dy(1)*(2*x(1) - 2*y(:,1)).^2)./(4*R0.^3) - dy(1)./R0 + (dy(2)*(2*x(1) - 2*y(:,1)).*(2*x(2) - 2*y(:,2)))./(4*R0.^3)) - dy(2)*((dy(2)*(2*x(2) - 2*y(:,2)).^2)./(4*R0.^3) - dy(2)./R0 + (dy(1)*(2*x(1) - 2*y(:,1)).*(2*x(2) - 2*y(:,2)))./(4*R0.^3));
%                         
%                     case 3
%                         error('Havent coded derivtives this high yet')
%                 end
%                 %alternative option if the sides meet at a vertex, or are
%                 %close together
%             end
%         end
%         
%         function R = distAnalCorner(self, sDist, t, t2corner, deriv, sSide, tSide)
%             % s basically means collocation point here, and
%             % t is the variable which moves through the support
%             theta = self.internalAngle(sSide,tSide);
%             Rdist = @(r) sqrt(sDist^2 + (t2corner + r).^2 - 2*cos(theta)*sDist*(t2corner + r));
%             
%             if sSide == tSide
%                 error('This method is designed for neighbouring sides')
%             else  
%                 R0 = Rdist(t);
%                 switch deriv
%                     case 0 
%                         R = R0;
%                     case 1
%                         R = (-sDist*cos(theta) + t2corner + t)./R0;
%                     case 2
%                         R = 1./R0 - (-2*sDist*cos(theta) + 2*(t2corner + t)).^2./(4*R0.^3);
%                     case 3
%                         error('Havent coded derivtives this high yet')
%                 end
%             end
%         end
    end
end

