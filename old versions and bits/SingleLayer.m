classdef SingleLayer < BIEoperator
    %single layer operator defined on a single domain
    
    properties
    end
    
    methods
        function self=SingleLayer(kwave,domain)
            self=self@ BIEoperator(kwave, @SingleLayer.SLkernel);
            %self.kernelInputs={'x','y','R'};
        end
    end
        
    methods (Static)
        function vals=SLkernel(kwave,x,y,varargin)
            if nargin==3 || isempty(varargin{1}) || isempty(varargin{1}{1}) %ad infinitum, ideally
                R=sqrt( (x(:,1)-y(:,1)).^2 + (x(:,2)-y(:,2)).^2);
            else
                R=varargin{1};
            end
            vals=(1i/4)*besselh(0,1,kwave*R);
        end
    end
    
    
end

