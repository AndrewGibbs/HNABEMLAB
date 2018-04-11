classdef singleLayer
    %single layer kernel
    
    %** so far everything assumes that x(s) and y(t) are on the same side
    
    properties
        kwave
        domain
        singularity = 'log'
        phase
%         phaseAnalExtDerivs
    end
    
    methods
        function self = singleLayer(kwave,domain)
%             if nargin == 2
%                 analDerivs=3;
%             end
            self.domain=domain;
            self.kwave=kwave;
            %now make the phase and its derivatives
            %sgn = @(s,t) sign(real(s-t));
            %differentiating w.r.t 2nd variable 't':
            self.phase = @(s,t) self.domain.dist(s,t);
%             for n=1:analDerivs
%                 self.phaseAnalExtDerivs{n}=@(s,t) domain.distAnal...
%             end
        end
        
        function g = phaseAnalDeriv(self,s,t,deriv,sGEt)
            if nargin == 5
                g = self.domain.distAnal(s,t,deriv,sGEt);
            else
                g = self.domain.distAnal(s,t,deriv);
            end
        end
        
        function K = kernel(self,s,t)
            %single layer kernel
            K = 1i/4 * besselh(0,1,self.kwave*abs(s-t));
        end
        
        function K = kernelNonOscAnal(self,s,t,sGEt)
            %non-oscillatory part of single layer kernel, extended
            %analytically (provided sgn is constant)
            if nargin == 4
                K = 1i/4 * besselh(0,1,self.kwave*self.domain.distAnal(s,t,0,sGEt))./exp(1i*self.kwave*self.domain.distAnal(s,t,0,sGEt));
            else
                K = 1i/4 * besselh(0,1,self.kwave*self.domain.distAnal(s,t,0))./exp(1i*self.kwave*self.domain.distAnal(s,t,0));
            end
        end
        
    end
end