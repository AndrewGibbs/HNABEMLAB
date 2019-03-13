classdef (Abstract) basis <handle
    %have children HNA and hp
    
    properties
        el
        obstacle
        mesh
        pMax
        nLayers
        sigmaGrad
        meshDOFs
        edgeBasis
        elSide
        plusCoefs = []
        minusCoefs = []
        nonOscCoeffs = []
    end
    
    methods
        function V=plus(V1,V2)
            %should introduce a check here for matching elements, to be
            %remvoed
            V.numEls=V1.numEls+V2.numEls;
            V.els=[V1.els V2.els];
            if isequal(V1.side,V2.side)
                V.side=V1.side;
            else
                %with this approach, understanding which basis element
                %points to which side is lost.
                V.side={V1.side,V2.side};
            end
        end
        
        %can't overload indexing function, as this already has a meaning
        function n=numEls(self)
            n=length(self.el);
        end
        
        function proj=project(self,coeffs)
            proj=Projection(coeffs,self);
        end
        
        function [proj,coeffs]=leastSquares(self,s,fs)
            %approximation at points s, of function
            %fs, by elements
            X=zeros(self.numEls,length(s));
            for j=1:self.numEls
                X(j,:)=self.el(j).eval(s);
            end         
            X=X.';
            coeffs=(X.'*X)\X.'*fs;
            proj=self.project(coeffs);
        end
    end
    
end