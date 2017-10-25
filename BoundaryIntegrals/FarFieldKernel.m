classdef FarFieldKernel < kernel
    %far field kernel
    
    properties
    end
    
    methods
        function  self = FarFieldKernel(kwave)
              self.kwave=kwave;
        end
        
        function [Val, OscVal, NonOscVal]=eval(self,theta,y)
            size_y=size(y);
            if size_y(2)~=2
                error('second input to far-field kernel must be an Nx2 vector');
            end
            theta=theta(:);
            NonOscVal=1;
            [Y1, CosTheta]=meshgrid(y(:,1), cos(theta));
            [Y2, SinTheta]=meshgrid(y(:,2), sin(theta));
            OscVal=exp(-1i*self.kwave*(Y1.*CosTheta+Y2.*SinTheta));
            Val=OscVal;
        end
    end
    
end

