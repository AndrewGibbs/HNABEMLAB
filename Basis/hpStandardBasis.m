classdef hpStandardBasis <basis
    
    properties
        %in parent basis
        hMax
    end
    
    methods
        function self=hpStandardBasis(side, pMax, hMax, nLayers, sigmaGrad)
            %store key parameters
            self.pMax=pMax;
            self.hMax=hMax;
            self.nLayers=nLayers;
            self.sigmaGrad=sigmaGrad;
            %first construct the mesh
            mesh=singleMesh(side, nLayers, sigmaGrad, hMax);
            
            self.meshDOFs=zeros(1,length(mesh1.el));
            
            %now put basis elements on it
            elCount=0;
            
            %initiate vector of objects, for some reason must be done this
            %way
            self.el=baseFnHNA(0,pMax,mesh.el(1),0,side);
                        
            for m=1:length(mesh.el)
                mesh.el(m).pMax=pMaxChoose( mesh.el(m).gradIndex, pMax, mesh.minDex );
                for p=0:mesh.el(m).pMax
                    elCount=elCount+1;
                    self.el(elCount)=baseFnHNA(0,p,mesh.el(m),0,side);
                    self.meshDOFs(m)=mesh.el(m).pMax+1;
                end
            end
            
%            self.numEls=elCount;
            self.side=side;
            self.mesh=mesh;
        end
    end
    
end

