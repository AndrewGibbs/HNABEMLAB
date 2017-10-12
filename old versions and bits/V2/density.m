classdef density < BoundaryFunction
    %general density function, designed for beam source problems
    
    properties
        evalFn
    end
    
    methods
        function self=density(domain, evalFn)
            self.domain=domain;
            self.supp=domain.supp;
            self.suppWidth=domain.L;
            if nargin==1
                self.evalFn=@(s) 1;
            else
                self.evalFn=evalFn;
            end
        end
        
        function val=eval(self,s)
           if min(s)<self.supp(1) || max(s)>self.supp(2)
               error('Input value outside of parametrisation range');
           else
               val=self.evalFn(s);
           end
        end
    end
    
end

