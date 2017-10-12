classdef HNAsingleMesh <basis
    
    properties
        %in parent basis
    end
    
    methods
        function self=HNAsingleMesh(side, pMax, kwave, alphaDist, nLayers, sigmaGrad)
            %first construct the mesh
            mesh=singleMesh(side, nLayers, sigmaGrad);
            
            %now put basis elements on it
            elCount=0;
            alphaDistUsed=0;
            
            %initiate vector of objects, for some reason must be done this
            %way
            self.el=baseFnHNA(kwave,pMax,mesh.el(1),0);
            
            for meshEl=mesh.el
                for p=0:pMaxChoose( meshEl.gradIndex, pMax, mesh.minDex )
                    if min(meshEl.distL, meshEl.distR)>=alphaDist*2*pi/kwave
                        for pm=[-1 1]
                            elCount=elCount+1;
                            self.el(elCount)=baseFnHNA(kwave,p,meshEl,pm);
                        end
                    else
                        alphaDistUsed=1;
                        elCount=elCount+1;
                        self.el(elCount)=baseFnHNA(kwave,p,meshEl,0);
                    end
                end
            end
            
            if alphaDistUsed==0
               warning('No basis elements removed, may lead to ill conditioning of discrete system '); 
            end
            self.numEls=elCount;
            self.side=side;
        end
    end
    
end

