classdef BEMintegral3D < BEMintegral
    %assumes two kernels, at least one of which is smooth. Splits into (potentially singular) 2D
    %integral and a smooth 1D integral
    %middle integral typically integrates over an entire obstacle
    
    properties
        singularKerIndex    %index of singular kernel
        Fn1
        Fn2
    end
    
    methods
        function self=BEMintegral3D(domain, kernel, Fn1, Fn2)
            self.Fn1=Fn1;
            self.Fn2=Fn2;
            self.domain=domain;
            self.kernel=kernel;
            
            if isequal(domain{1},domain{2}) && isequal(domain{2},domain{3})
                error('Both kernels may be singular, code cannot handle this yet');
            end
            if isequal(domain{1},domain{2})
                self.singularKerIndex=1;
            elseif isequal(domain{2},domain{3})
                self.singularKerIndex=2;
            end    
            self.type='tricky...';
            
        end
    end
    
end