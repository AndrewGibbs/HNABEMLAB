classdef HNAsingleMesh <basis
    
    properties
        %in parent basis
        alphaDist
        plusCoefs = []
        minusCoefs = []
        nonOscCoeffs = []
        edgeBasis
        elSide
    end
    
    methods
        function self=HNAsingleMesh(obstacle, pMax, kwave, alphaDist, nLayers, sigmaGrad, oscNearEndPoints)
            %store key parameters
            if nargin<=6
                oscNearEndPoints=0;
            end
                
            self.pMax=pMax;
            self.nLayers=nLayers;
            self.sigmaGrad=sigmaGrad;
            self.alphaDist=alphaDist;
            self.obstacle=obstacle;
            
            if isa(obstacle,'polygon')
                indexEnd = 0;
                self.meshDOFs = [];
                for n=1:obstacle.numSides
                    indexStart = indexEnd + 1;
                    self.edgeBasis{n} = HNAsingleMesh(obstacle.side{n}, pMax, kwave, alphaDist, nLayers, sigmaGrad, oscNearEndPoints);
                    indexEnd = indexStart -1 + length(self.edgeBasis{n}.el);
                    self.mesh{n} = self.edgeBasis{n}.mesh;
                    
                    self.plusCoefs    = [self.plusCoefs    indexStart+self.edgeBasis{n}.plusCoefs];
                    self.minusCoefs   = [self.minusCoefs   indexStart+self.edgeBasis{n}.minusCoefs];
                    self.nonOscCoeffs = [self.nonOscCoeffs indexStart+self.edgeBasis{n}.nonOscCoeffs];
                    self.el           = [self.el           self.edgeBasis{n}.el];
                    self.meshDOFs     = [self.meshDOFs     self.edgeBasis{n}.meshDOFs];
                    self.elSide(indexStart:indexEnd) = n;
                end
                return;
            end
            
            %first construct the mesh
            mesh=singleMesh(obstacle, nLayers, sigmaGrad);
            
            self.meshDOFs=zeros(1,length(mesh.el));
            
            %now put basis elements on it
            elCount=0;
            alphaDistUsed=0;
            
            %initiate vector of objects, for some reason must be done this
            %way
            self.el=baseFnHNA(kwave,pMax,mesh.el(1),0,obstacle);
            
            for m=1:length(mesh.el)
                mesh.el(m).pMax=pMaxChoose( mesh.el(m).gradIndex, pMax, mesh.minDex );
                
                %first establish which type of phase we want, given
                %distance to the corners, and the alpha parameter for how
                %much should be thrown away:
                if mesh.el(m).width>=alphaDist*2*pi/kwave 
                    PM = [-1 1];
                else
                    alphaDistUsed=1;
                    if oscNearEndPoints==1
                        if mesh.el(m).distL<mesh.el(m).distR %closer to left hand endpoint
                            PM=1;
                        else
                            PM=-1;
                        end
                    else
                        PM=0;
                    end
                end
                %previously compared against: min(mesh.el(m).distL,
                %mesh.el(m).distR). But have adjusted because a)
                %this is kind of wrong, and b) to agree with
                %Emile's thesis.
                for pm=PM
                    for p=0:mesh.el(m).pMax
                        elCount=elCount+1;
                        self.el(elCount)=baseFnHNA(kwave,p,mesh.el(m),pm,mesh.side);
                        mesh.el(m).osc=1;
                        self.meshDOFs(m)=self.meshDOFs(m)+1;
                    end
                end
            end
            
            
            %now store which coeffs are which. Initialise the coefficient vectors
            self.plusCoefs=[];
            self.minusCoefs=[];
            self.nonOscCoeffs=[];
            for n=1:length(self.el)
                switch self.el(n).pm
                    case 1
                        self.plusCoefs = [self.plusCoefs n];
                    case -1
                        self.minusCoefs = [self.minusCoefs n];
                    case 0
                        self.nonOscCoeffs = [self.nonOscCoeffs n];
                end
            end
            
            if alphaDistUsed==0
%               warning('No basis elements removed, may lead to ill conditioning of discrete system UNLESS singular values are removed '); 
            end
%            self.numEls=elCount;
            self.mesh=mesh;
            
            self.elSide = ones(size(self.el));
            
        end
    end
   
    
end

