classdef (Abstract) mesh
    
    properties
        points
        el
        side
        minDex
    end
    
    methods
        
        function flipMesh=uminus(self)
            %maybe need copy here:
            flipMesh=self;
            flipMesh.points=sort(max(self.points)-self.points);
            for m=1:self.numEls;
                flipMesh.el(m).interval=[flipMesh.points(m)  flipMesh.points(m+1)];
                flipMesh.el(m).gradIndex=self.el(self.numEls-m+1).gradIndex;
                flipMesh.el(m).width=self.el(self.numEls-m+1).width;
            end
            
             %now construct corner distances
             for m=1:length(self.el)
                 self.el(m).distL=0;
                 for m_=1:(m-1) %have to manually sum
                    self.el(m).distL=self.el(m).distL+sum(self.el(m_).width);
                 end
                self.el(length(self.el)-m+1).distR=sum(self.el(m).distL);
             end
        end
        
        function hMax=maxWidth(self)
            hMax=0;
            for n=1:length(self.el)
                if self.el(n).width>hMax
                    hMax=self.el(n).width;
                end
            end
        end
        
       function n=numEls(self)
            n=length(self.points) -1;
        end
    end
    
end