classdef (Abstract) BoundaryFunction < handle
    properties
       domain
       edgeComponent
       numEdgeComponents
       suppEdges
    end
    
    methods 
        
        function self = BoundaryFunction(domain, edgeComponents, suppEdges)
           if nargin == 2 
               self.suppEdges = true( domain.numComponents,1);
           else
               self.suppEdges = suppEdges;
           end
           self.numEdgeComponents = domain.numComponents;
           self.domain = domain;
           for n=1:self.numEdgeComponents
               self.edgeComponent{n} = edgeComponents(self.domain.component{n});
           end
        end
        
        function val = eval(self, x, nEdge)
            val = self.edgeComponent{nEdge}.eval(x);
        end
        
        function val = nonOscAnal(self, x, nEdge)
            val = self.edgeComponent{nEdge}.nonOscAnal(x);
        end
        
        function val = phaseAnal(self, x, deriv, nEdge)
            val = self.edgeComponent{nEdge}.phaseAnal(x,deriv);
        end
    end
    
end


%         function ResBasFn=restrictTo(self,ResDomain,ResWidth)
% %             ResBasFn=copy(self);
% %             ResBasFn.supp=ResDomain;
%           ResBasFn=self;
%             if nargin==2
%                 ResBasFn.suppWidth=ResDomain(2)-ResDomain(1);
%             else
%                 ResBasFn.suppWidth=ResWidth;
%             end
%           ResBasFn.supp=ResDomain;
%         end
%         
%         function I=L2(f,g)
%             if isa(f,'BoundaryIntegral') && length(f.domain)==2 && isa(g,'BoundaryFunction')
%                 I=BEMintegral2D(f.domain, f.kernel, g, f.boundaryFn);
%             elseif isa(f,'BoundaryFunction')
%                 I=InnerProduct1D( f, g);
%             end
%         end
%         
%         function supp = getSupp(self,side)
%             if isa(self.domain,'edge')
%                 supp = self.supp;
%             elseif isa(self.domain,'polygon')
%                 supp = self.supp(side,:);
%             end
%         end