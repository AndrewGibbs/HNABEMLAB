classdef (Abstract) BoundaryFunction 
    %IF THIS EVER BECOMES A HANDLE CLASS, NEED TO CHANGE RESTRICTION BELOW
    %have children HNA and hp
    
    % ** could add 'a' and 'b' properties here, for convinience
    
    % ** should ultimately include a nonOscAnal bit in here, to force all
    % subclasses to share this method.
    properties
        supp %mesh interval over which fn is supported
        suppWidth %measure of support
        domain %side on which the support lives, when parametrised, supp is a subset
        oscillator %required for any oscillatory integration
        %phase
        a
        b
        phaseMaxStationaryPointOrder
    end
    
    methods 
        eval(obj)
        nonOscAnal(obj)
        phaseAnal(obj)
        
        function ResBasFn=restrictTo(self,ResDomain,ResWidth)
%             ResBasFn=copy(self);
%             ResBasFn.supp=ResDomain;
          ResBasFn=self;
            if nargin==2
                ResBasFn.suppWidth=ResDomain(2)-ResDomain(1);
            else
                ResBasFn.suppWidth=ResWidth;
            end
          ResBasFn.supp=ResDomain;
        end
        
        function I=L2(f,g)
            if isa(f,'BoundaryIntegral') && length(f.domain)==2 && isa(g,'BoundaryFunction')
                I=BEMintegral2D(f.domain, f.kernel, g, f.boundaryFn);
            elseif isa(f,'BoundaryFunction')
                I=InnerProduct1D( f, g);
            end
        end
        
        function supp = getSupp(self,side)
            if isa(self.domain,'edge')
                supp = self.supp;
            elseif isa(self.domain,'polygon')
                supp = self.supp(side,:);
            end
        end
    end
    
end