classdef (Abstract) waveR2
    %incident field
    
    properties
        kwave
        phaseLinear
    end
    
    methods
        val=DirTrace(self,s,boundary);
        val=NeuTrace(self,s,boundary);
        val=POA(self,s,boundary);
        val=phase(self,s,boundary);
    end
    
end