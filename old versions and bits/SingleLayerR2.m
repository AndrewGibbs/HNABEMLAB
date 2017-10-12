classdef SingleLayerR2 < waveR2 & BoundaryIntegral
    %single layer operator defined on R^2
    
    properties
        %domain
        %Fn %essentially the density in this case
    end
    
    methods
        function self=SingleLayerR2(kwave,domain,density,integratorDef)
            if nargin==3
                self.integratorDef=[];
            else
                self.integratorDef=integratorDef;
            end
            self.kwave=kwave;
            self.domain=domain;
            self.boundaryFn=density;
            self.kernel=@(R) (1i/4)*besselh(0,1,self.kwave*R);
            self.NeuKernel=@(R,n_xms) NeuKernel(self.kwave,R,n_xms);
            %for domain and disjoint integrals
%            self.type='smooth';
            self.toR2=1;
            self.integratorDef;
            %now create altrnate versions of kernels needed for integration
        end
    end
        
    
    methods
        function val=eval(self,x,integrator)
            xSize=size(x);
            if xSize(2)~=2
                error('first input must be Nx2 vector');
            end
            if nargin==2
                if isequal(self.integrator,[])
                    error('Integrator must be specified');
                else
                    val=self.integratorDef.eval(self,x);
                end
            else
                val=integrator.eval(self,x);
            end
        end
        
        function val=DirTrace(self, boundary, s, integrator)
            %need to check that domains are disjoint here
            if isequal(boundary,self.domain)
                error('Can only take trace to surface away from source');
            end
            s=s(:);
            S=SingleLayer(self.kwave,{self.domain, boundary});
            Sg=S*self.density;
            if nargin==3
                if isequal(self.integrator,[])
                    error('Integrator must be specified');
                else
                    val=self.integratorDef.eval(Sg,x);
                end
            else
                val=integrator.eval(Sg,x);
            end
            
        end
        
        function val=NeuTrace(self,boundary,s)
            %need to check that domains are disjoint here
            I=boundaryIntegral(self.NeuKernel,self.boundaryFn,boundary);
            val=integrator(I,s);
        end
       
    end
    
end

