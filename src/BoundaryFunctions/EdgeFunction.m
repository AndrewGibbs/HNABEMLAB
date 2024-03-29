classdef (Abstract) EdgeFunction 
    %IF THIS EVER BECOMES A HANDLE CLASS, NEED TO CHANGE RESTRICTION BELOW
    %have children HNA and hp
    
    % ** should ultimately include a nonOscAnal bit in here, to force all
    % subclasses to share this method.
    properties
        supp %mesh interval over which fn is supported
        suppWidth %measure of support
        domain %side on which the support lives, when parametrised, supp is a subset
        oscillator %required for any oscillatory integration
        a
        b
        phaseMaxStationaryPointOrder
        meshEl
        nodesPerWavelength = 20 %was 20
    end
    
    methods 
        eval(obj)
        nonOscAnal(obj)
        phaseAnal(obj)
        
        function S = getSupp(self,side)
            if isa(self.domain,'edge')
               S = self.supp; 
            end
        end
    end
    
end