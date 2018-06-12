classdef hpStandardBasis <basis
    
    properties
        %in parent basis
        hMax
        edgeBasis
        elSide
    end
    
    methods
        function self=hpStandardBasis(obstacle, pMax, hMax, nLayers, sigmaGrad)
            %store key parameters
            self.pMax=pMax;
            self.hMax=hMax;
            self.nLayers=nLayers;
            self.sigmaGrad=sigmaGrad;
            self.obstacle=obstacle;
            
            
            if isa(obstacle,'polygon')
                indexEnd = 0;
                self.meshDOFs = [];
                for n=1:obstacle.numSides
                    indexStart = indexEnd + 1;
                    self.edgeBasis{n} = hpStandardBasis(obstacle.side{n}, pMax, hMax, nLayers, sigmaGrad);
                    indexEnd = indexStart -1 + length(self.edgeBasis{n}.el);
                    self.mesh{n} = self.edgeBasis{n}.mesh;
                    
%                     self.plusCoefs    = [self.plusCoefs    indexStart+self.edgeBasis{n}.plusCoefs];
%                     self.minusCoefs   = [self.minusCoefs   indexStart+self.edgeBasis{n}.minusCoefs];
%                     self.nonOscCoeffs = [self.nonOscCoeffs indexStart+self.edgeBasis{n}.nonOscCoeffs];
                    self.el           = [self.el           self.edgeBasis{n}.el];
                    self.meshDOFs     = [self.meshDOFs     self.edgeBasis{n}.meshDOFs];
                    self.elSide(indexStart:indexEnd) = n;
                end
                return;
            end
            
            %first construct the mesh
            mesh=singleMesh(obstacle, nLayers, sigmaGrad, hMax);
            
            self.meshDOFs=zeros(1,length(mesh.el));
            
            %now put basis elements on it
            elCount=0;
            
            %initiate vector of objects, for some reason must be done this
            %way
            self.el=baseFnHNA(0,pMax,mesh.el(1),0,obstacle);
                        
            for m=1:length(mesh.el)
                mesh.el(m).pMax=pMaxChoose( mesh.el(m).gradIndex, pMax, mesh.minDex );
                for p=0:mesh.el(m).pMax
                    elCount=elCount+1;
                    self.el(elCount)=baseFnHNA(0,p,mesh.el(m),0,obstacle);
                    self.meshDOFs(m)=mesh.el(m).pMax+1;
                end
            end
            
%            self.numEls=elCount;
            self.mesh=mesh;
            self.elSide = ones(size(self.el));
        end
    end
    
end