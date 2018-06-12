classdef HNAoverlappingMesh <basis
    
    properties
        %in parent basis
    end
    
    methods
        function self=HNAoverlappingMesh(obstacle, pMax, kwave, nLayers, sigmaGrad)
            %store key parameters
            self.pMax=pMax;
            self.nLayers=nLayers;
            self.sigmaGrad=sigmaGrad;
            
            %first construct the mesh
            mesh1=meshSide(obstacle, nLayers, sigmaGrad);
            mesh2=-mesh1;
            
            self.meshDOFs{1}=zeros(1,length(mesh1.el));
            self.meshDOFs{1}=zeros(1,length(mesh2.el));
            
            self.el=baseFnHNA(kwave,pMax,mesh1.el(1),0,mesh1.side);
            elCount=0;
            
            for m=1:length(mesh1.el)
                mesh1.el(m).pMax=pMaxChoose( mesh1.el(m).gradIndex, pMax, mesh1.minDex );
                for p=0:mesh1.el(m).pMax
                    elCount=elCount+1;
                    self.el(elCount)=baseFnHNA(kwave,p,mesh1.el(m),1,mesh1.side);
                    self.meshDOFs{1}(m)=mesh1.el(m).pMax+1;
                end
            end
            for m=1:length(mesh2.el)
                mesh2.el(m).pMax=pMaxChoose( mesh2.el(m).gradIndex, pMax, mesh2.minDex );
                for p=0:mesh2.el(m).pMax
                    elCount=elCount+1;
                    self.el(elCount)=baseFnHNA(kwave,p,mesh2.el(m),-1,mesh2.side);
                    self.meshDOFs{2}(m)=mesh2.el(m).pMax+1;
                end
            end
            
%            self.numEls=elCount;
            self.obstacle=obstacle;
            self.mesh{1}=mesh1;
            self.mesh{2}=mesh2;
        end
    end
    
end