classdef HNAoverlappingMesh <basis
    
    properties
        %in parent basis
    end
    
    methods
        function self=HNAoverlappingMesh(side,pMax,kwave, nLayers, sigmaGrad)
            %first construct the mesh
            mesh1=meshSide(side, nLayers, sigmaGrad);
            mesh2=-mesh1;
            
            self.el=baseFnHNA(kwave,pMax,mesh1.el(1),0,mesh1.side);
            elCount=0;
            
            for meshEl=mesh1.el
                for p=0:pMaxChoose( meshEl.gradIndex, pMax, mesh1.minDex )
                    elCount=elCount+1;
                    self.el(elCount)=baseFnHNA(kwave,p,meshEl,1,mesh1.side);
                end
            end
            for meshEl=mesh2.el
                for p=0:pMaxChoose( meshEl.gradIndex, pMax, mesh1.minDex )
                    elCount=elCount+1;
                    self.el(elCount)=baseFnHNA(kwave,p,meshEl,-1,mesh1.side);
                end
            end
            
%            self.numEls=elCount;
            self.side=side;
            self.mesh{1}=mesh1;
            self.mesh{2}=mesh2;
        end
    end
    
end