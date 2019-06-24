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
            self.domain=domain;
            self.kwave=kwave;
        end
        
        function g = phaseAnalDeriv(self,s,t,deriv,sGEt,sSide,tSide)
            g = self.domain.distAnal(s,t,deriv,sGEt,sSide,tSide);
        end
        
        function g = phaseAnalDerivCorner( self, sDist, t, t2corner, deriv, sSide, tSide)
            g = self.domain.distAnalCorner( sDist, t, t2corner, deriv, sSide, tSide);
        end
        
        function K = kernel(self, s, t, sSide, tSide)
            %single layer kernel
            if nargin == 2
                %radial kernel, where s = |x-y|
                K = 1i/4 * besselh(0,1,self.kwave*s);
                return;
            end
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
        
        function K = kernelNonOscAnalCorner(self, sDist, t, t2corner, sSide, tSide)
            %non-oscillatory part of single layer kernel, extended for
            %corner cases
%                 K = 1i/4 * besselh(0,1,self.kwave*self.domain.distAnalCorner( sDist, t, t2corner, 0, sSide, tSide))...
%                     ./exp(1i*self.kwave*self.domain.distAnalCorner( sDist, t, t2corner, 0, sSide, tSide));       
                [~,~,K_] = besselhDecomp(0,1,self.kwave*self.domain.distAnalCorner( sDist, t, t2corner, 0, sSide, tSide));
                K = K_*1i/4;
        end
        
        function m = phaseMaxStationaryPointOrder(self,sameSide)
            % ** need to add some extra conditions in here when I do
            % polygons
            if isa(self.domain,'Screen') || isa(self.domain,'MultiScreen')
               m=0; 
            elseif isa(self.domain,'ConvexPolygon')
                if ~sameSide
                    m=1;
                else
                    m=0;
                end
            else
               error('Domain type not recognised');
            end
        end
        function X = getSymmetries(self)
            X = abs(self.domain.getSymmetries());
        end
    end
end