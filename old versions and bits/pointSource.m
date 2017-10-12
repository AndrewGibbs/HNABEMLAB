classdef pointSource < incField
    %plane wave incident field
    
    properties
        s
    end
    
    methods
        
        function self=pointSource(kwave,direction)
           self.kwave=kwave;
           self.d=direction;
        end
        
        function val=dirTrace(self,t,boundary)
            %standard def of point source
            val=(1i/4)*besselh(0,1i*k*self.s*boundary.trace(t));
        end
        
        function val=neuTrace(self,t,boundary)
            %taken from old code
            x=poly.trace(t); n_x=boundary.n_s(t);
            R=sqrt((poly.source(1)-x(:,1)).^2+(poly.source(2)-x(:,2)).^2);
            n_dot_x_minus_source=(n_x(:,1).*x(:,1) + n_x(:,2).*x(:,2))...
                                  - (poly.source(1).*n_x(:,1) + poly.source(2).*n_x(:,2));
            val=-(1i*k/4)*besselh(1,1,self.kwave*R).*n_dot_x_minus_source./R;
        end
        
        function val=POA(self,t,boundary)
            %taken from old code
            n_x=boundary.n_s(t);
            P_x= boundary.vertex(boundary.arc_side(t),:); %get the corners in terms of t
            illum_s=((n_x(:,1).*P_x(:,1) + n_x(:,2).*P_x(:,2))-(poly.n_s(t)*poly.source'))<0;
            val=2*neuTrace(self,t,boundary).*illum_s;
        end
    end
    
end

