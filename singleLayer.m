classdef singleLayer
    %single layer kernel
    
    %** so far everything assumes that x(s) and y(t) are on the same side
    
    properties
        kwave
        obstacle
        singularity = 'log'
        phase
    end
    
    methods
        function self = singleLayer(kwave,obstacle)
            self.obstacle=obstacle;
            self.kwave=kwave;
            %now make the phase and its derivatives
            sgn = @(s,t) sign(real(s-t));
            self.phase = {@(s,t) sgn(s,t)*(s-t), @(s,t) sgn(s,t), @(s,t) 0};
        end
        
        function K = kernel(s,t)
            %single layer kernel
            K = 1i/4 * besselh(0,1,self.kwave*abs(s-t));
        end
        
        function K = kernelNonOscAnal(s,t)
            %non-oscillatory part of single layer kernel, extended
            %analytically (provided sgn is constant)
            sgn = sign(real(s-t));
            K = 1i/4 * besselh(0,1,self.kwave*sgn*(s-t))./exp(1i*self.kwave*sgn*(s-t));
        end
        
    end
end

