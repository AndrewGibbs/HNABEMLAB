classdef ProjectionEdge < EdgeFunction
    
    properties
        coeffs
        edgeBasis
        sideDOFs
        nodes
        weights
    end
    
    methods
        function self = ProjectionEdge(component,coeffs,basis)
            %find the index/side of the component of 'component'
            for n=1:length(basis.obstacle.component)
               if isequal(basis.obstacle.component(n),component)
                   compIndex = n;
                   break;
               end
            end
            self.edgeBasis = basis.edgeBasis{n};
            self.coeffs = coeffs(basis.elEdge==compIndex);
            self.sideDOFs = length(self.coeffs);
            self.domain = component;
            
           k = self.edgeBasis.el(1).kwave;
           [self.nodes, self.weights] = self.edgeBasis.getPoints(self.nodesPerWavelength,k,'G');
        end
        
        function val = eval(self,s)
            
           s=s(:);
           val=zeros(length(s),1);
           for j=1:self.sideDOFs
               val=val+self.edgeBasis.el(j).eval(s)*self.coeffs(j);
           end
        end
        
    end
end

