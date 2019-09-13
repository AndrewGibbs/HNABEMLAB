classdef HNAoverlappingMesh <basis
    
    properties
        %in parent basis
        plusCoefs =[]
        minusCoefs =[]
    end
    
    methods
        function self=HNAoverlappingMesh(obstacle, pMax, kwave, nLayers, sigmaGrad)
            %store key parameters
            self.pMax=pMax;
            self.nLayers=nLayers;
            self.sigmaGrad=sigmaGrad;
            self.obstacle=obstacle;
            
            if ~isa(obstacle,'edge')
                MultiBasis(@(x) HNAoverlappingMesh(x, pMax, kwave, nLayers, sigmaGrad), obstacle, self);
                return;
            end
            
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
                    self.plusCoefs = [self.plusCoefs elCount];
                end
            end
            for m=1:length(mesh2.el)
                mesh2.el(m).pMax=pMaxChoose( mesh2.el(m).gradIndex, pMax, mesh2.minDex );
                for p=0:mesh2.el(m).pMax
                    elCount=elCount+1;
                    self.el(elCount)=baseFnHNA(kwave,p,mesh2.el(m),-1,mesh2.side);
                    self.meshDOFs{2}(m)=mesh2.el(m).pMax+1;
                    self.minusCoefs = [self.minusCoefs elCount];
                end
            end
            
            self.mesh{1}=mesh1;
            self.mesh{2}=mesh2;
            
%             if isa(obstacle,'polygon')
%                 error('Overlapping mesh not yet compatible with polygons, use single mesh');
%             else
%                 self.edgeBasis = []; %not sure what this does yet...
%                 self.elSide(1:length(self.el))=1;
%             end
        end
        
        function [s, w] = getPoints(self,pointsPerWavelength,k,distType)
            if self.sigmaGrad>=.5
                error('Need smaller value of sigmaGrad to obtain mesh info');
            end
            %returns points with fixed number of nodes per mesh element,
           fakeSingleMesh = self.mimicSingleMesh;
           [s,w] = getQuadPointsFromMesh(fakeSingleMesh,pointsPerWavelength,k,distType);
            
        end
        
        function [fakeMesh, DOFs] = mimicSingleMesh(self)
            %merge the two overlapping meshes to make a single mesh, which
            %is useful sometimes, for allocating collocation & quadrature
            fakeMesh.points = [self.mesh{1}.points(1:(end-1)) self.mesh{2}.points(2:end)];%unique(sort([self.mesh{1}.points self.mesh{2}.points]));
            %fakeMesh.side = component;
            %numMeshEls = length(meshPoints)-1;
            %get widths
            for n = 1:(length(self.mesh{1}.el)-1)
                  fakeMesh.el(n) = self.mesh{1}.el(n);
                  DOFs(n) = self.meshDOFs{1}(n);
            end
            
            %now merge the two middle mesh elements:
            halfMeshWidthsLength = length(fakeMesh.el);
            
            fakeMesh.el(halfMeshWidthsLength+1).width = abs(self.mesh{2}.el(1).interval(2) - self.mesh{1}.el(end).interval(1));
            fakeMesh.el(halfMeshWidthsLength+1).interval = [self.mesh{1}.el(end).interval(1) self.mesh{2}.el(1).interval(2)];
            fakeMesh.el(halfMeshWidthsLength+1).distL = self.mesh{1}.el(halfMeshWidthsLength+1).distL;
            fakeMesh.el(halfMeshWidthsLength+1).distR = fakeMesh.el(halfMeshWidthsLength+1).distL;
            %self.mesh{1}.el(n).distR;
            
            fakeMesh.el(halfMeshWidthsLength+1).gradIndex = self.mesh{1}.el(n).gradIndex;
            fakeMesh.el(halfMeshWidthsLength+1).pMax = self.mesh{1}.el(n).pMax;
            DOFs(halfMeshWidthsLength+1) = self.meshDOFs{1}(end) + self.meshDOFs{2}(1);
            
            for n = 2:(length(self.mesh{2}.el))
                fakeMesh.el(halfMeshWidthsLength+n) = self.mesh{2}.el(n);
                DOFs(halfMeshWidthsLength+n) = self.meshDOFs{2}(n);
            end
        end
    end
    
end