classdef (Abstract) kernel
    %kernel class, designed to be used in conjunction with BoundaryOperator
    
    properties 
        kwave
        coeff=1
    end
    
    methods        
        
        function kernelTimes=mtimes(c,self)
            kernelTimes=self;
            kernelTimes.coeff=c;
        end
    end
    
end