classdef (Abstract) incField
    %incident field
    
    properties
        kwave
    end
    
    methods
        val=dirTrace(self,boundary);
        val=neuTrace(self,boundary)
    end
    
end

