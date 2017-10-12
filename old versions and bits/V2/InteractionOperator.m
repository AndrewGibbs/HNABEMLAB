classdef InteractionOperator < BoundaryOperator
    
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function self=InteractionOperator(fromDom, toDom)
            self.domain{1}=fromDom;
            self.domain{2}=toDom;
        end
        
        function val=kernel(s,t)
            %for mappings to anything other than a screen, need to 
           
            R=sqrt((Z1-Y1).^2 + (Z2-Y2).^2);
            val=-0.5*(1i*k*(besselh(1,1,k*R).*NdYmZ.*R)./(R.^2));
        end
    end
    
end

