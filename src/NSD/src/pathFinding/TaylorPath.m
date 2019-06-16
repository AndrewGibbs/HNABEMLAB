classdef TaylorPath
    %approximates h_{\xi,m}(p) using Taylor series
    
    properties
        coeffs
    end
    
    methods
        function self = TaylorPath(SPorder, G, startPoint, branch, ascFlag)
            if nargin <=4
                ascFlag = false;
            end
            %use botched version of the IC finder to get Taylor coeffs:
            h0Derivs = NSDpathICvAllTerms( SPorder, G, startPoint, ascFlag );
            %get Taylor coeffs of the higher order terms, these are just
            %the initial conditions from the ODE
            for j = 1:length(h0Derivs{branch})
                self.coeffs(j) = h0Derivs{branch}(j);
            end
        end
        
        function hp = path(self,p,deriv)
            if nargin == 2
                deriv = 0;
            end
            
            hp = 0;
            for n=0:(length(self.coeffs) - deriv - 1)
                hp = hp + self.coeffs(n+1+deriv)*p.^n/factorial(n);
            end
        end
    end
end

