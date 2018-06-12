classdef pointSource < waveR2
    %plane wave incident field
    
    properties
        source
    end
    
    methods
        
        function self=pointSource(kwave,sourcePoint)
           self.kwave=kwave;
           self.source=sourcePoint;
           self.phaseMaxStationaryPointOrder=1;
        end
        
        function [valOsc, valNonOsc]=DirTrace(self,t,boundary)
            %standard def of point source
            x=boundary.trace(t);
            R=sqrt((self.source(1)-x(:,1)).^2+(self.source(2)-x(:,2)).^2);
            valOsc=(1i/4)*besselh(0,1,self.kwave*R);
            valNonOsc=valOsc.*exp(-1i*self.kwave*R);
        end
        
        function [valOsc, valNonOsc]=NeuTrace(self,t,boundary)
            %taken from old code
            x=boundary.trace(t); n_x=boundary.normal(t);
            R=sqrt((self.source(1)-x(:,1)).^2+(self.source(2)-x(:,2)).^2);
            n_dot_x_minus_source=(n_x(:,1).*x(:,1) + n_x(:,2).*x(:,2))...
                                  - (self.source(1).*n_x(:,1) + self.source(2).*n_x(:,2));
            valOsc=-(1i*self.kwave/4)*besselh(1,1,self.kwave*R).*n_dot_x_minus_source./R;
            valNonOsc=valOsc.*exp(-1i*self.kwave*R);
        end
        
        function g = phasePD(self,X,xDers,yDers)
            dist = sqrt((X(:,1)-self.source(1)).^2 + (X(:,2)-self.source(2)).^2);
            %now list partial derivatives of the distance function
            if xDers == 0 && yDers == 0
                g = dist;
            elseif xDers == 1 && yDers == 0
                g = (X(:,1)-self.source(1))./dist;
            elseif xDers == 0 && yDers == 1
                g = (X(:,2)-self.source(2))./dist;
            elseif xDers == 1 && yDers == 1
                g = -(X(:,1)-self.source(1)).*(X(:,2)-self.source(2))./(dist.^3);
            elseif xDers == 2 && yDers == 0
                g = (X(:,1)-self.source(1)).^2./(dist.^3);
            elseif xDers == 0 && yDers == 2
                g = (X(:,2)-self.source(2)).^2./(dist.^3);
            else
                error('Havent coded such high derivatives for distance function');
            end
        end
        
        function x = sourceVsnormal(self,side)
            x = -sign((side.P1-self.source)*(side.nv.'));
        end
    end
    
end