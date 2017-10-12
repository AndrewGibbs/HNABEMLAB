classdef SLkernel < kernel
    %kernel class, designed to be used in conjunction with BoundaryOperator
    
    properties 
        inputs={'r'}
    end
    
    methods        
        
        function self=SLkernel(kwave)
            self.kwave=kwave;
        end
        
        function val=eval(self,r)
           val=self.coeff*(1i/4)*besselh(0,1,self.kwave*r); 
        end
        
        function val=evalFull(self,s,t,S,T,matFlag)
           if nargin==5
               matFlag=0;
           end
           x=S.trace(s); y=T.trace(t);
           if matFlag==0
                R=sqrt((x(:,1)-y(:,1)).^2 + (x(:,2)-y(:,2)).^2);
           else
               [X1, Y1]=meshgrid(x(:,1),y(:,1));
               [X2, Y2]=meshgrid(x(:,2),y(:,2));
               R=sqrt((X1-Y1).^2+(X2-Y2).^2);
           end
           val=self.eval(R);
        end
    end
    
end