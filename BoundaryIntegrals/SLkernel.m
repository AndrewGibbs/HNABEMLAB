classdef SLkernel < kernel
    %kernel class, designed to be used in conjunction with BoundaryOperator
    
    properties 
        inputs={'r'}
    end
    
    methods        
        
        function self=SLkernel(kwave)
            self.kwave=kwave;
        end
        
        function [val, nonOscVal] =eval(self,r)
           val=self.coeff*(1i/4)*besselh(0,1,self.kwave*r); 
           nonOscVal=val./exp(1i*self.kwave*r);
        end
        
        function  [val, nonOscVal]=evalFull(self,s,t,S,T,matFlag)
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
           [val, nonOscVal]=self.eval(R);
        end
        
    end
    
    methods(Static)
        
        function linFun=colPhaseLin(x,xGLy)
            if nargin==1
                linFun=[];
            elseif xGLy=='xLy'
                linFun=[1,-x];
            elseif xGLy=='yLx'
                linFun=[-1,x];
            else
                error('Unrecognised format xGLy in kernel phase');
            end
        end
    end
    
end