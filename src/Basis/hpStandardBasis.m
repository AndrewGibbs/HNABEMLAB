classdef hpStandardBasis <basis
    
    properties
        %in parent basis
        hMax
    end
    
    methods
        function self=hpStandardBasis(obstacle, pMax, hMax, nLayers, sigmaGrad)
            %store key parameters
            self.pMax=pMax;
            self.hMax=hMax;
            self.nLayers=nLayers;
            self.sigmaGrad=sigmaGrad;
            self.obstacle=obstacle;
            
            if ~isa(obstacle,'edge')
                MultiBasis(@(x) hpStandardBasis(x, pMax, hMax, nLayers, sigmaGrad), obstacle, self);
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
            
            self.plusCoefs = [];
            self.minusCoefs = [];
%            self.numEls=elCount;
            self.nonOscCoeffs=1:elCount;
            self.mesh=mesh;
%            self.elEdge = ones(size(self.el));
        end
        
        function [s, w] = getPoints(self,pointsPerWavelength,k,distType)
            if self.sigmaGrad>=.5
                error('Need smaller value of sigmaGrad to obtain mesh info');
            end
            %returns points with fixed number of nodes per mesh element,
           fakeSingleMesh = self.mesh;
           k = ceil(pi/self.hMax);
           [s,w] = getQuadPointsFromMesh(fakeSingleMesh,pointsPerWavelength,k,distType);
            
        end
    end
    
end