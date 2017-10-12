classdef DLkernel < kernel
    %kernel class, designed to be used in conjunction with BoundaryOperator
    
    properties 
        inputs={'r','n_dot_x_minus_y'}
    end
    
    methods        
        
        function self=DLkernel(kwave)
            self.kwave=kwave;
        end
        
        function val=eval(self,r,n_dot_x_minus_y)
            val=-self.coeff*(1i*self.kwave/4)*besselh(1,1,self.kwave*r).*n_dot_x_minus_y./r;
        end
        
        function val=evalFull(self,s,t,S,T,matFlag)
            if nargin==5
                matFlag=0;
            end
           x=S.trace(s); y=T.trace(t); nx=S.nv;
           %[X1,Y1]=meshgrid(x(:,1),y(:,1));[X2,Y2]=meshgrid(x(:,2),y(:,2));
           if matFlag==0
                R=sqrt((x(:,1)-y(:,1)).^2 + (x(:,2)-y(:,2)).^2);
                NdXmY=nx(1)*(x(:,1)-y(:,1)) + nx(2)*(x(:,2)-y(:,2));
           else
               [X1, Y1]=meshgrid(x(:,1),y(:,1));
               [X2, Y2]=meshgrid(x(:,2),y(:,2));
               R=sqrt((X1-Y1).^2+(X2-Y2).^2);
               NdXmY=nx(1)*(X1-Y1) + nx(2)*(X2-Y2);
           end
           val=self.eval(R,NdXmY);
        end
    end
    
end