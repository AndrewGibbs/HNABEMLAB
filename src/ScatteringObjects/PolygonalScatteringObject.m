classdef (Abstract) PolygonalScatteringObject < handle
    %abstract class of scattering objects, such as screens, convex
    %polygons, etc
    
    properties
        component
        numComponents
        vertices
        internalAngle
        Lipschitz
    end
    
    methods
        val = normal(self,s);
        
        function draw(self)
           plot(self.vertices(:,1),self.vertices(:,2)) ;
        end
    
        function val = trace(self,s,index)
            val = self.component{index}.trace(s);
        end
        
        function val = dist(self,s,sIndex,t,tIndex)
            val = abs(self.trace(s,sIndex)-self.trace(t,tIndex));
        end
         
        function R = distAnal(self,s,t,deriv,sIndex,tIndex)
            if length(s)>1
                error('first input must be a single param value');
            end
            s=s(:); t=t(:);
            if sIndex == tIndex
                x = self.side{sIndex}.trace(s);
                y = self.side{tIndex}.trace(t);
                yMx = -[ (x(1)-y(:,1))  (x(2)-y(:,2)) ];
                dy = self.side{tIndex}.dSv;
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
                        R = - dy(1)*((dy(1)*(2*x(1) - 2*y(:,1)).^2)./(4*R0.^3) - dy(1)./R0 + (dy(2)*(2*x(1) - 2*y(:,1)).*(2*x(2) - 2*y(:,2)))./(4*R0.^3)) - dy(2)*((dy(2)*(2*x(2) - 2*y(:,2)).^2)./(4*R0.^3) - dy(2)./R0 + (dy(1)*(2*x(1) - 2*y(:,1)).*(2*x(2) - 2*y(:,2)))./(4*R0.^3));
                        
                    case 3
                        error('Havent coded derivtives this high yet')
                end
                %alternative option if the sides meet at a vertex, or are
                %close together
            end
        end
    end
end

