classdef (Abstract) waveR2
    %incident field
    
    properties
        kwave
        oscillator  %oscillator of derivative should be the same
    end
    
    methods
        val=DirTrace(self,s,boundary);
        val=NeuTrace(self,s,boundary);
        vsl=POA(self,s,boundary);
    end
    
end