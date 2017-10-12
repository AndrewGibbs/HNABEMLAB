classdef BoundaryIntegral < BoundaryFunction
    %general class which deals with 'Af(x)' type operator-function
    %relations.
    %can also describe 'A...Bf(x)' type operator-functions
    
    properties
        domain
        kernel
        Fn
        toR2          %flag, 0 if  maps to another boundary, 1 if to R2?
        integratorDef=[] %default integration routine
    end
    
    methods 
        function self=BoundaryIntegral(kernel,boundaryFn,boundary)
            if nargin==2
                self.toR2=1;
                boundary=[];
            else
                self.toR2=0;
            end
            
            if isa(boundaryFn,'BoundaryIntegral');
                self.kernel=horzcat(boundaryFn.kernel,kernel);
                self.domain=horzcat(boundaryFn.domain,boundary);
                self.Fn=boundaryFn.boundaryFn;
            else
                self.kernel=kernel;
                self.domain=boundary;
                self.Fn=boundaryFn;
            end
        end
        
        function out=mtimes(m,self)
            out=BoundaryIntegral(m*self.kernel,self.Fn,self.domain);
        end
        
        function setIntegrator(self,integratorIn)
            self.integratorDef=integratorIn;
        end
%         function val=ptEval(self,pt,integrator)
%             I=BEMintegral1D(self.kernel, funciton, pt);
%             val=I.eval(integrator);
%         end
        
        function val=eval(self,X)
            if isequal(self.integratorDef,[])
                error('Integrator needs to be set');
            else
                val=self.integratorDef(self,X);
            end
        end
        
%         function val=innerProduct(self,boundaryFn2,integrator)
%             if length(self.domain)==2
%                 I=BEMintegral2D(self.kernel, self.Fn, boundaryFn2);
%                 val=I.eval(integrator);
%             elseif length(self.domain)==3
%                 I=BEMintegral3D(self.kernel, self.Fn, boundaryFn2);
%                 val=I.eval(integrator);
%             else
%                 error('Cant go to higher dimensional integrals than this');
%             end
%         end
        
        function I=L2(Kf,g)
            if ~(isa(Kf,'BoundaryIntegral') && isa(g,'BoundaryFunction'))
                error('First arguemnt must be BoundaryIntegral, second BoundaryFunction');
            end
            if length(Kf.domain)==2
                I=BEMintegral2D(Kf.domain,Kf.kernel, Kf.Fn, g);
            elseif length(Kf.domain)==3
                I=BEMintegral3D(Kf.domain, Kf.kernel, Kf.Fn, g);
            end
        end
    end
    
end