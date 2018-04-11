classdef pointSource < waveR2
    %plane wave incident field
    
    properties
        source
    end
    
    methods
        
        function self=pointSource(kwave,sourcePoint)
           self.kwave=kwave;
           self.source=sourcePoint;
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
        
%         function val=POA(self,t,boundary)
%             %taken from old code
%             n_x=boundary.n_s(t);
%             P_x= boundary.vertex(boundary.arc_side(t),:); %get the corners in terms of t
%             illum_s=((n_x(:,1).*P_x(:,1) + n_x(:,2).*P_x(:,2))-(poly.n_s(t)*poly.source'))<0;
%             val=2*neuTrace(self,t,boundary).*illum_s;
%         end
    end
    
end

