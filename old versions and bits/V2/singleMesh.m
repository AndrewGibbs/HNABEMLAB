classdef singleMesh  < mesh
    
    properties
       %all inhereted 
    end
    
    methods
        
        function self=singleMesh(side, nLayers, sigmaGrad, hMax)
            if nargin<=2
                sigmaGrad=0.15;
            end
            if nargin==3
                hMax=inf;
            end
            if sigmaGrad>=.5
                error('grading param must be less than .5');
            end
            self.side=side;
            
            %create two meshes without splitting for hMax
            meshL=meshSide(side, nLayers, sigmaGrad);
            meshR=-meshL;
            self.minDex=meshL.minDex;
            %create copy as starting point for 'self'
            %self=meshL;
            
            self.points=sort(union(meshL.points,meshR.points));
            
            %create new middle element
            newMidEl.interval=[meshL.points(end-1) meshR.points(2)];
            newMidEl.gradIndex=meshL.el(end).gradIndex;
            newMidEl.width=newMidEl.interval(2)-newMidEl.interval(1); %relatively large, so this move is OK
            
            %define these bits as NaNs for now, so we can concatenate
            newMidEl.distL=NaN;
            newMidEl.distR=NaN;
            
            
            %put it all together
            self.el=[meshL.el(1:(end-1)) newMidEl meshR.el(2:end)];
            
            %now split elements based on hMax
            m=1;
            while hMax<self.maxWidth
                 %elPreSplit=self.el;
                 %for m=1:length(self.el)
                    if hMax<inf && self.el(m).width>hMax
                        splits=ceil( self.el(m).width/hMax); %split interval this many times
                        for n=1:splits
                            splitEls(n).interval=[self.el(m).interval(1)+(n-1)*self.el(m).width/splits self.el(m).interval(1)+n*self.el(m).width/splits   ];
                            splitEls(n).gradIndex=self.el(m).gradIndex;
                            splitEls(n).width=self.el(m).width/splits;
                            %and just bodge this bit, will get fixed shortly
                            splitEls(n).distL=NaN;
                            splitEls(n).distR=NaN;
                        end
                        self.el=[self.el(1:(m-1)) splitEls self.el((m+1):end)];
                        clear splitEls;
                    end
                 %end
                 m=m+1;
            end
            
%            self.numEls=length(self.el);
            
             %now construct corner distances
             for m=1:length(self.el)
                 self.el(m).distL=0;
                 for m_=1:(m-1) %have to manually sum
                    self.el(m).distL=self.el(m).distL+sum(self.el(m_).width);
                 end
                self.el(length(self.el)-m+1).distR=sum(self.el(m).distL);
             end
                
            
            %reconstruct 'points' of mesh given new elements
            self.points=zeros(1,self.numEls+1);
            for m=1:self.numEls
                self.points(m)=self.el(m).interval(1);
            end
            self.points(end)=self.el(end).interval(2);            

        end
        
    end
    
end

