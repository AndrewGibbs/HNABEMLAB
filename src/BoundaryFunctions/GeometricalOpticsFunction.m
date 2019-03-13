classdef GeometricalOpticsFunction < BoundaryFunction
    %only valid for screens!
    
    properties
        uinc
    end
    
    methods
        function self=GeometricalOpticsFunction(uinc,domain)
            %self.edgeComponent = GeometricalOpticsEdge(uinc,edge());
            
            %+-1 to determine direction of reflection
            if ~domain.Lipschitz %non-Lipschitz domain
                dirConst(1:domain.numComponents) = -uinc.sourceVsnormal(domain);
                illumEdges = 1:domain.numComponents;
            else %on a polygon, zero on sides which can't 'see' inc wave
                illumEdges = [];
                for n=1:domain.numComponents
                    pre_dirConst = uinc.sourceVsnormal(domain.component(n));
                    if pre_dirConst<0
                        dirConst(n) = 1;
                        illumEdges = [illumEdges n];
                    else
                        dirConst(n) = 0;
                    end
                    clear pre_dirConst;
                end
            end
            
            edgeComponents = @(x) GeometricalOpticsEdge(uinc,x);
            
            self@BoundaryFunction(domain, edgeComponents, illumEdges);
            
            self.uinc=uinc;
            self.domain=domain;
            
            for n=1:domain.numComponents
                self.edgeComponent(n).dirConst = dirConst(n);
            end
        end
    end
    
end