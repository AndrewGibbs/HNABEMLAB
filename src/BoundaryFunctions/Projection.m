classdef Projection < BoundaryFunction
    %projection onto a discrete space, via coefficients and a basis
    
    properties
        coeffs
        plusCoefs
        minusCoefs
        el
        elSide
        Ndim
        GOA
        onScreen
        nodes
        weights
    end
    
    methods
        function self=Projection(coeffs,basis)
            
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
           self.plusCoefs = basis.plusCoefs;
           self.minusCoefs = basis.minusCoefs;
           if isa(basis.obstacle,'edge')
               self.onScreen=true;
           else
               self.onScreen=false;
           end
           
           k = basis.el(1).kwave;
           [self.nodes, self.weights] = basis.mesh.getPoints(self.nodesPerWavelength,k,'G');
        end  
        
        F = FarField( self, theta );
        
        function val=eval(self,s, side)
            if nargin==2
                if self.onScreen
                    side=1;
                else
                    error('Need to specify which side please.');
                end
            end
            s=s(:);
            val=zeros(length(s),1);
            allDofs = 1:self.Ndim;
            sideDofs = allDofs(self.elSide == side);
           for j=sideDofs
               val=val+self.el(j).eval(s)*self.coeffs(j);
           end
        end
        
        function val=evalComp(self, s, pm, side)
            if nargin==3
                if self.onScreen
                    side=1;
                else
                    error('Need to specify which side please.');
                end
            end
            s=s(:);
            val=zeros(length(s),1);
            allDofs = 1:self.Ndim;
            sideDofs = allDofs(self.elSide == side);
            if pm=='+'
                subCoeffs = self.plusCoefs;
            elseif pm=='-'
                subCoeffs = self.minusCoefs;
            else
                error('third arg needs to be char + or -');
            end
            sideSubDofs = intersect(sideDofs,subCoeffs);
            val=zeros(length(s),1);
            
           for j=sideSubDofs
               val=val+self.el(j).eval(s)*self.coeffs(j);
           end
        end
        
        function v_pm = getDiffEnv(self,pm)
                v_pm = DiffApprox(self,pm);
        end
    end
    
end