classdef (Abstract) baseFn 
    %IF THIS EVER BECOMES A HANDLE CLASS, NEED TO CHANGE RESTRICTION BELOW
    %have children HNA and hp
    properties
        supp %mesh interval over which basis fn is supported
        side %side on which the support lives
    end
    
    methods 
        val(obj)
        
        function ResBasFn=restrictTo(self,ResDomain)
%             ResBasFn=copy(self);
%             ResBasFn.supp=ResDomain;
              ResBasFn=self;
              ResBasFn.supp=ResDomain;
        end
    end
    
end