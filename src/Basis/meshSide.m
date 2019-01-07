%not sure if the mesh splitting option is needed here at all.

classdef meshSide < mesh
    %general mesh class
    
    properties
        %
        %sigmaGrad
        %maxIndex
        %preSplitMeshPoints
    end
    
    methods
        function self=meshSide(side, nLayers, sigmaGrad, hMax)
            if nargin<=2
                sigmaGrad=0.15;
            end
            if nargin==3
                hMax=inf;
            end
            
            if mod(nLayers,1)~=0 || nLayers<0
                error('Number of mesh layeres must be non-negative integer');
            end
            %self.sigmaGrad=sigmaGrad;
            normalisedPoints=sort([0 sigmaGrad.^(0:nLayers)]);
            self.points=side.L*normalisedPoints;
            %self.preSplitMeshPoints=self.points;
            meshIndices=1:length(self.points) -1;
            %self.maxIndex=max(meshIndices);
            
            widths1=side.L*sort([sigmaGrad^nLayers sigmaGrad.^((1:nLayers)-1)*(1-sigmaGrad)]);
            %split each element such that width<hMax
            for m=meshIndices
                if hMax<inf && widths1(m)>hMax
                    splits=ceil( widths1(m)/hMax); %split interval this many times
                    widths=[widths1(1:(m-1))  ones(1,splits)*widths1(m)/splits  widths1((m+1):end)];
                    self.points=[self.points(1:m) (self.points(m)+(1:(splits-1))*widths(m)) self.points((m+1):end)];
                    meshIndices=[meshIndices(1:(m-1)) ones(1,splits)*meshIndices(m) meshIndices((m+1):end)];
                else
                    widths=widths1;
                end
            end
            %self.numEls=length(self.points) -1;
            
            for m=[self.numEls 1:(self.numEls-1)]
                self.el(m).interval=[self.points(m)  self.points(m+1)];
                self.el(m).gradIndex=meshIndices(m); %usful later for grading polynomial degree
                self.el(m).width=widths(m);
            end
            self.side=side;
            
             %now construct corner distances
             for m=1:length(self.el)
                 self.el(m).distL=0;
                 for m_=1:(m-1) %have to manually sum
                    self.el(m).distL=self.el(m).distL+sum(self.el(m_).width);
                 end
                %self.el(length(self.el)-m+1).distR=sum(self.el(m).distL);
                self.el(m).distR=side.L-self.el(m).distL-self.el(m).width;
             end
             
             %now define the 'minDex' variable for graded polynomial degree
             self.minDex=1;
             while self.el(m).interval(1)<2 && self.minDex<self.numEls
                 self.minDex=self.minDex+1;
             end
        end
    end
    
end