classdef BoundaryOperator < handle
    %general BIE class. Not quite abstract.
    
    properties
        kernel
        kwave
        domain
    end
    
    methods
        function self=BIEoperator(kwave,kernel,domain)
             if nargin==3
                 self.codomain=domain;
             else 
                 self.codomain=codomain;
             end
             self.domain=domain;
             self.kernel=@(kwave,x,y,varargin) kernel(kwave,x,y,varargin);
             self.kwave=kwave;
        end

        %parametrised kernel function
        function kernel_st(self, side_s, s, side_t, t, R)
            if nargin==6
                self.kernel([],[],R);
            elseif isequal(side_s,side_t)
                R=abs(s-t);
                self.kernel([],[],R);
            else
                x=side_s.trace(s); y=side_t.trace(t);
                self.kernel(x,y);
            end
        end
            
        
        function K=plus(A,B)
            if isequal(A.kwave, B.kwave)
                kernelSum=@(x,y,varargin) A.kernel(x,y,varargin) + B.kernel(x,y,varargin);
                K=BIEoperator(A.kwave,kernelSum);
            else
                error('wavenumber of operators dont match, cant be added.');
            end
        end
        
        function K=mtimes(c,A)
            %need a check that 'c' is 1D here
            if isequal(size(c),[1,1])
                kernelScaled=@(kwave,x,y,extras) c*A.kernel(kwave,x,y,extras);
                K=BIEoperator(A.kwave,kernelScaled);
            else
                error('cannot multiple operator by vector');
            end
            
        end
    end
    
end

