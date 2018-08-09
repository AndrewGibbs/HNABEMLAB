classdef polygon < scatteringObject
    %contains all information and methods needed about a polygon
    
    properties
        side
        cumL
        extensionJoinL
        extensionJoinR
        internalAngle
    end
    
    methods
        function self = polygon(vertices)
            [sizeY, sizeX]= size(vertices);
            %do some basic sense checks
%             if sizeY < 3
%                 error('A polygon must have at least 3 vertices');
%             end
            if sizeX ~= 2
                error('Input must be Nx2 vector');
            end
            %add first vertex on the end for good measure:
            self.numSides = sizeY;
            self.vertices = [vertices; vertices(1,:)];
            
            self.L=zeros(self.numSides,1);
            
            %create each side of the polygon as an 'edge'
            for j = 1:sizeY
                self.side{j} = edge(self.vertices(j:(j+1),:));
                self.L(j) = self.side{j}.L;
                self.side{j}.nv = - self.side{j}.nv;
            end
            self.supp = [0 sum(self.L)];
            self.cumL=[0; cumsum(self.L)];
            
            %now make matrix of internal angles,
            %and distances of how far the side must be extended to meet
            %another side
            for n = 1:self.numSides
                for m= 1: self.numSides
                    if n == m
                        self.internalAngle(m,n) = pi;
                        extensionJoin = [0; 0];
                    elseif norm(self.side{m}.dSv - self.side{n}.dSv )<10^-12
                        %sides are parallel
                        self.internalAngle(m,n) = NaN;
                        extensionJoin = [NaN; NaN];
                    else
                        extensionJoin =  intersect_lines( self.side{n}.dSv.', self.side{n}.P1.', self.side{m}.dSv.', self.side{m}.P1.' );%...
%                         [self.side{n}.dSv(1) self.side{m}.dSv(1); 
%                          self.side{n}.dSv(2) self.side{m}.dSv(2)]...
%                          \(self.side{n}.P2.' - self.side{m}.P1.');
                        self.internalAngle(m,n) = acos(- self.side{m}.dSv * (self.side{n}.dSv.'));
                    end
                    self.internalAngle(n,m) = self.internalAngle(m,n);
                    self.extensionJoinL(n,m) = extensionJoin(1);
                    self.extensionJoinR(n,m) = extensionJoin(2);
                end
            end
            
        end
        
        function [sideOf_s, s_onSide ] = getSide(self,s)
            %given the parametrisation of the polygon, what side are we on?
            sideOf_s = NaN(size(s));
            for j = 1:self.numSides
               onSide = self.cumL(j) <= s & s < self.cumL(j+1);
               s_onSide{j} = s(onSide);
               sideOf_s(onSide) = j;
            end
        end
        
        function val = trace(self,s,ofSide)
            val = self.side{ofSide}.trace(s);
        end
        
        function R = distAnal(self,s,t,deriv,sGEt,sSide,tSide)
            s=s(:); t=t(:);
            if sSide == tSide
                R = self.side{sSide}.distAnal(s,t,deriv,sGEt);
            else  % some condition about the sides not being too close together...
                x = self.side{sSide}.trace(s);
                y = self.side{tSide}.trace(t);
                yMx = -[ (x(1)-y(:,1))  (x(2)-y(:,2)) ];
                dy = self.side{tSide}.dSv;
                R0 = sqrt( (x(1) - y(:,1) ).^2 + (x(2) - y(:,2) ).^2);
                switch deriv
                    case 0 
                        R = R0;
                    case 1
                        R = (dy*(yMx.')).'./R0;
                    case 2
                        %account for 2nd derivative bug:
                        if length(t) == 1
                            yMx1 = yMx(1);
                            yMx2 = yMx(2);
                        else
                            yMx1 = yMx(1,:);
                            yMx2 = yMx(2,:);
                        end
                      %  R = (dy * ( [dy(1)./R0 dy(2)./R0] - [y(:,1)./(R0.^3)  y(:,2)./(R0.^3)] * sum(yMx1*dy(1) + yMx2*dy(2)) ).').';
                        R = - dy(1)*((dy(1)*(2*x(1) - 2*y(:,1)).^2)./(4*R0.^3) - dy(1)./R0 + (dy(2)*(2*x(1) - 2*y(:,1)).*(2*x(2) - 2*y(:,2)))./(4*R0.^3)) - dy(2)*((dy(2)*(2*x(2) - 2*y(:,2)).^2)./(4*R0.^3) - dy(2)./R0 + (dy(1)*(2*x(1) - 2*y(:,1)).*(2*x(2) - 2*y(:,2)))./(4*R0.^3));
                        
                    case 3
                        error('Havent coded derivtives this high yet')
                end
                %alternative option if the sides meet at a vertex, or are
                %close together
            end
        end
        
        function R = distAnalCorner(self, sDist, t, t2corner, deriv, sSide, tSide)
            % s basically means collocation point here, and
            % t is the variable which moves through the support
            theta = self.internalAngle(sSide,tSide);
            Rdist = @(r) sqrt(sDist^2 + (t2corner + r).^2 - 2*cos(theta)*sDist*(t2corner + r));
            
            if sSide == tSide
                error('This method is designed for neighbouring sides')
            else  
                R0 = Rdist(t);
                switch deriv
                    case 0 
                        R = R0;
                    case 1
                        R = (-sDist*cos(theta) + t2corner + t)./R0;
                    case 2
                        R = 1./R0 - (-2*sDist*cos(theta) + 2*(t2corner + t)).^2./(4*R0.^3);
                    case 3
                        error('Havent coded derivtives this high yet')
                end
            end
        end
    end
end

