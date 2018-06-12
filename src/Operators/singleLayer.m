classdef singleLayer
    %single layer operator, contains everything you need to know about it
    
    %** so far everything assumes that x(s) and y(t) are on the same side
    
    properties
        kwave
        domain
        singularity = 'log'
        %phase
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
            %self.phase = @(s,t) self.domain.dist(s,t);
%             for n=1:analDerivs
%                 self.phaseAnalExtDerivs{n}=@(s,t) domain.distAnal...
%             end
        end
        
        function g = phaseAnalDeriv(self,s,t,deriv,sGEt,sSide,tSide)
%             if nargin == 5
                g = self.domain.distAnal(s,t,deriv,sGEt,sSide,tSide);
%             else
%                 g = self.domain.distAnal(s,t,deriv);
%             end
        end
        
        function K = kernel(self,s,t, sSide, tSide)
            %single layer kernel
            if  sSide == tSide
                K = 1i/4 * besselh(0,1,self.kwave*abs(s-t));
            else
                K = 1i/4 * besselh(0,1,self.kwave*abs(self.domain.distAnal(s,t,0,[],sSide,tSide)));
            end
        end
        
        function K = kernelNonOscAnal(self,s,t,sGEt, sSide, tSide)
            %non-oscillatory part of single layer kernel, extended
            %analytically (provided sgn is constant)
            if nargin == 6
                K = 1i/4 * besselh(0,1,self.kwave*self.domain.distAnal(s,t,0,sGEt, sSide, tSide))./exp(1i*self.kwave*self.domain.distAnal(s,t,0,sGEt, sSide, tSide));
            else
                K = 1i/4 * besselh(0,1,self.kwave*self.domain.distAnal(s,t,0))./exp(1i*self.kwave*self.domain.distAnal(s,t,0));
            end
            
        end
        
        function m = phaseMaxStationaryPointOrder(self,sameSide)
            % ** need to add some extra conditions in here when I do
            % polygons
            if isa(self.domain,'edge')
               m=0; 
            elseif isa(self.domain,'polygon')
                if ~sameSide
                    m=1;
                else
                    m=0;
                end
            else
               error('Domain type not recognised');
            end
        end
    end
end