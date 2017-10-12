classdef hpStandardBasis <basis
    
    properties
        %in parent basis
    end
    
    methods
        function self=hpStandardBasis(side, pMax, hMax, nLayers, sigmaGrad)
            %first construct the mesh
            mesh=singleMesh(side, nLayers, sigmaGrad, hMax);
            
            %now put basis elements on it
            elCount=0;
            
            %initiate vector of objects, for some reason must be done this
            %way
            self.el=baseFnHNA(0,pMax,mesh.el(1),0,side);
            
            for meshEl=mesh.el
                for p=0:pMaxChoose( meshEl.gradIndex, pMax, mesh.minDex )
                    elCount=elCount+1;
                    self.el(elCount)=baseFnHNA(0,p,meshEl,0,side);
                end
            end
            
%            self.numEls=elCount;
            self.side=side;
            self.mesh=mesh;
        end
    end
    
end

