classdef InnerProduct1D < BEMintegral
    %inner product, no kernel
    
    properties
        Fn1
        Fn2
        supp
    end
    
    methods
        function self=InnerProduct1D(Fn1,Fn2)
           self.Fn1=Fn1; 
           self.Fn2=Fn2; 
           self.kernel=[];
           if isequal(Fn1.domain,Fn2.domain)
               self.supp=intersect( Fn1.supp, Fn2.supp );
               self.domain=Fn1.domain;
           end
        end
    end
    
end