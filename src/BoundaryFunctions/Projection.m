classdef Projection < BoundaryFunction
    %projection onto a discrete space, via coefficients and a basis
    
    properties
        coeffs
        el
        elSide
        Ndim
        GOA
    end
    
    methods
        function self=Projection(coeffs,basis, GOA)
            if nargin==2
                %add an option GOA to the solution
                self.GOA=[];
            else
                self.GOA=GOA;
            end
            if basis.numEls~=length(coeffs)
                error('Dimensions of basis and coefficients vector do not match');
            end
            self.coeffs=coeffs;
           self.el=basis.el; 
           self.elSide=basis.elSide; 
           self.domain=basis.obstacle;
           self.supp=basis.obstacle.supp;
           self.suppWidth=basis.obstacle.L;
           self.Ndim=basis.numEls;
        end  
        
        F = FarField( self, theta );
        
        function val=eval(self,s, side)
            s=s(:);
            val=zeros(length(s),1);
            allDofs = 1:self.Ndim;
            sideDofs = allDofs(self.elSide == side);
           for j=sideDofs
               val=val+self.el(j).eval(s)*self.coeffs(j);
           end
        end
    end
    
    
end