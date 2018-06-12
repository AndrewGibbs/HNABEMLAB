classdef (Abstract) waveR2
    %incident field
    
    properties
        kwave
        phaseLinear
        phaseMaxStationaryPointOrder
    end
    
    methods
        val=DirTrace(self,s,boundary);
        val=NeuTrace(self,s,boundary);
        val=POA(self,s,boundary);
        val=phasePD(self,s,xDers,yDers);
        val=sourceVsnormal(self,edge);
        %val=phase(self,s,boundary);
    end
    
end