classdef Projection < BoundaryFunction
    %projection onto a discrete space, via coefficients and a basis
    
    properties
        coeffs
        el
        Ndim
    end
    
    methods
        function self=Projection(coeffs,basis)
            if basis.numEls~=length(coeffs)
                error('Dimensions of basis and coefficients vector do not match');
            end
            self.coeffs=coeffs;
           self.el=basis.el; 
           self.domain=basis.side;
           self.supp=[0 basis.side.L];
           self.suppWidth=basis.side.L;
           self.Ndim=basis.numEls;
        end  
        
        function val=eval(self,s)
            s=s(:);
            val=zeros(length(s),1);
           for j=1:self.Ndim
               val=val+self.el(j).eval(s)*self.coeffs(j);
           end
        end
    end
    
    
end