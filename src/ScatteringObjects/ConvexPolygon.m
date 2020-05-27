classdef ConvexPolygon < PolygonalScatteringObject
    %contains all information and methods needed about a polygon
    
    properties
        L %some measure of total length of boundary
    end
    
    methods
        function self = ConvexPolygon(vertices)
            
            [sizeY, sizeX]= size(vertices);
            if sizeX ~= 2
                error('Input must be Nx2 vector');
            end
            %add first vertex on the end for good measure:
            self.numComponents = sizeY;
            self.vertices = [vertices; vertices(1,:)];
            
            self.L=zeros(self.numComponents,1);
            
            %create each side of the polygon as an 'edge'
            self.L = 0;
            for j = 1:sizeY
                self.component(j) = edge(self.vertices(j:(j+1),:));
                self.L = self.L + self.component(j).L;
                self.component(j).nv = - self.component(j).nv;
            end
            
            %now make matrix of internal angles,
            %and distances of how far the side must be extended to meet
            %another side
            for n = 1:self.numComponents
                for m= 1: self.numComponents
                    if n == m
                        self.internalAngle(m,n) = pi;
                    elseif norm(self.component(m).dSv - self.component(n).dSv )<10^-12
                        %sides are parallel
                        self.internalAngle(m,n) = NaN;
                    else
                        self.internalAngle(m,n) = acos(- self.component(m).dSv * (self.component(n).dSv.'));
                    end
                    self.internalAngle(n,m) = self.internalAngle(m,n);
                end
            end
            
            self.Lipschitz = true;
            
        end
        
        function val = trace(self,s,ofSide)
            val = self.component(ofSide).trace(s);
        end
        
%         function R = distAnal(self,s,t,deriv,sGEt,sSide,tSide)
%             s=s(:); t=t(:);
%             if sSide == tSide
%                 R = self.component{sSide}.distAnal(s,t,deriv,sGEt);
%             else  % some condition about the sides not being too close together...
%                 x = self.component{sSide}.trace(s);
%                 y = self.component{tSide}.trace(t);
%                 yMx = -[ (x(1)-y(:,1))  (x(2)-y(:,2)) ];
%                 dy = self.component{tSide}.dSv;
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

