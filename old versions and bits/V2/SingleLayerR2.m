classdef SingleLayerR2 < waveR2 & BoundaryIntegral
    %single layer operator defined on R^2
    
    properties
%         domain
%         boundaryFn
        %Fn %essentially the density in this case
        isIntegral=1
    end
    
    methods
        function self=SingleLayerR2(kwave,densityFn,boundary)
            domain{1}='R2'; domain{2}=boundary;
            self@BoundaryIntegral(SLkernel(kwave),densityFn,domain);
            %self.kernel=@(R) (1i/4)*besselh(0,1,self.kwave*R);
            self.kwave=kwave;
            self.domain=domain;
            self.Fn=densityFn;
%            self.NeuKernel=@(R,n_xms) NeuKernel(self.kwave,R,n_xms);
            %for domain and disjoint integrals
%            self.type='smooth';
            self.toR2=1;
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
        
        function I=DirTrace(self, boundary)
            newDomain{1}=boundary;
            newDomain{2}=self.domain{2};
            I=BoundaryIntegral(self.kernel,self.Fn,newDomain);
        end
        
        function I=NeuTrace(self,boundary)
            newDomain{1}=boundary;
            newDomain{2}=self.domain{2};
            I=BoundaryIntegral(DLkernel(self.kwave),self.Fn,newDomain);
        end
       
    end
    
end