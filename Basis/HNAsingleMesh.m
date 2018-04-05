classdef HNAsingleMesh <basis
    
    properties
        %in parent basis
        alphaDist
    end
    
    methods
        function self=HNAsingleMesh(side, pMax, kwave, alphaDist, nLayers, sigmaGrad, oscNearEndPoints)
            %store key parameters
            if nargin<=6
                oscNearEndPoints=0;
            end
                
            self.pMax=pMax;
            self.nLayers=nLayers;
            self.sigmaGrad=sigmaGrad;
            self.alphaDist=alphaDist;
            
            
            %first construct the mesh
            mesh=singleMesh(side, nLayers, sigmaGrad);
            
            self.meshDOFs=zeros(1,length(mesh.el));
            
            %now put basis elements on it
            elCount=0;
            alphaDistUsed=0;
            
            %initiate vector of objects, for some reason must be done this
            %way
            self.el=baseFnHNA(kwave,pMax,mesh.el(1),0,side);
            
            for m=1:length(mesh.el)
                mesh.el(m).pMax=pMaxChoose( mesh.el(m).gradIndex, pMax, mesh.minDex );
                for p=0:mesh.el(m).pMax
                    if mesh.el(m).width>=alphaDist*2*pi/kwave 
                        %previously compared against: min(mesh.el(m).distL,
                        %mesh.el(m).distR). But have adjusted because a)
                        %this is kind of wrong, and b) to agree with
                        %Emile's thesis.
                        for pm=[-1 1]
                            elCount=elCount+1;
                            self.el(elCount)=baseFnHNA(kwave,p,mesh.el(m),pm,mesh.side);
                            self.meshDOFs(m)=self.meshDOFs(m)+1;
                        end
                        mesh.el(m).osc=1;
                    else
                        alphaDistUsed=1;
                        elCount=elCount+1;
                        if oscNearEndPoints==1
                            if mesh.el(m).distL<mesh.el(m).distR %closer to left hand endpoint
                                pm=1;
                            else
                                pm=-1;
                            end
                        else
                            pm=0;
                        end
                        self.el(elCount)=baseFnHNA(kwave,p,mesh.el(m),pm,mesh.side);
                        mesh.el(m).osc=0;
                        self.meshDOFs(m)=self.meshDOFs(m)+1;
                    end
                end
            end
            
            if alphaDistUsed==0
               warning('No basis elements removed, may lead to ill conditioning of discrete system '); 
            end
%            self.numEls=elCount;
            self.side=side;
            self.mesh=mesh;
        end
    end
    
end

