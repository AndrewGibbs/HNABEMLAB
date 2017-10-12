classdef PseudoBoundaryFunction < BoundaryFunction
    %takes a boundaryIntegral of the form Ag(x), and constructs an object
    %which can be treated as a boundaryFunction f=Ag
    
    properties
        integrator
        BoundaryIntegral
    end
    
    methods
        function self=PseudoBoundaryFunction(Ag,integrator)
            self.integrator=integrator;
            %BoundaryFunction type properties
            self.domain=Ag.domain{1};
            self.supp=self.domain.supp;
            self.suppWidth=self.domain.L;
            self.BoundaryIntegral=Ag;
        end
        
        function val=eval(self,X)
            val=self.integrator.eval(self.BoundaryIntegral,X);
        end
    end
    
end

