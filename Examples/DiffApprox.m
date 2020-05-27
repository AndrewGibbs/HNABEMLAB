classdef DiffApprox < BoundaryFunction
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        v_N %full diffracted 
        pm % {'+','-'}
        nodes
        weights
    end
    
    methods
        function self = DiffApprox(v_N, pm)
            self.nodes = v_N.nodes;
            self.weights = v_N.weights;
            self.v_N = v_N;
            self.pm = pm;
        end
        
        function vals = eval(self, s, side)
            
            if nargin==2
                if self.v_N.onScreen
                    side=1;
                else
                    error('Need to specify which side please.');
                end
            end
            
            vals = self.v_N.evalComp(s, self.pm, side);
        end
    end
end

