classdef combinedLayer < forumulation
    %combined layer operator, contains everything you need to know about it
    
    properties
%         kwave
%         domain
        coupling_param
        singularity = 'log'
    end
    
    methods
        function self = combinedLayer(kwave,domain,coupling_param)
            if nargin == 2
                self.coupling_param = kwave;
            else
                self.coupling_param = coupling_param;
            end
            self.domain=domain;
            self.kwave=kwave;
        end
        
        function g = phaseAnalDeriv(self,s,t,deriv,sGEt,sSide,tSide)
            g = self.domain.distAnal(s,t,deriv,sGEt,sSide,tSide);
        end
        
        function K = kernel(self, s, t, sSide, tSide)
            eta = self.coupling_param;
            %single layer kernel
            if nargin == 2
                %radial kernel, where s = |x-y|
                K = eta/4 * besselh(0,1,self.kwave*s);
                return;
            end
            if  sSide == tSide
                K = eta/4 * besselh(0,1,self.kwave*abs(s-t));
            else
                R = abs(self.domain.dist(s,t,sSide,tSide));
                n_x = self.domain.component(sSide).nv;
                n_dot_x_minus_y = (self.domain.component(sSide).trace(s) - self.domain.component(tSide).trace(t))*n_x.';
                K = -1i/4*self.kwave*(besselh(1,1,self.kwave*R).*n_dot_x_minus_y)./R + eta/4 * besselh(0,1,self.kwave*R);
            end
        end
        
        function K = kernelNonOscAnal(self,s,t,sGEt, sSide, tSide)
            %non-oscillatory part of single layer kernel, extended
            %analytically (provided sgn is constant)
            if sSide ~= tSide
                error('havent coded analytic version of this kernel for differnet sides yet');
            end
            if nargin == 6
                K = self.coupling_param/4 * besselh_0_1_nonosc_large_imag(self.kwave*self.domain.distAnal(s,t,0,sGEt, sSide, tSide));
                %self.coupling_param/4 * besselh(0,1,self.kwave*self.domain.distAnal(s,t,0,sGEt, sSide, tSide))./exp(1i*self.kwave*self.domain.distAnal(s,t,0,sGEt, sSide, tSide));
            else
                K = self.coupling_param/4 * besselh(0,1,self.kwave*self.domain.distAnal(s,t,0))./exp(1i*self.kwave*self.domain.distAnal(s,t,0));
            end
        end
        
        %identity component of operator:
        function Ig_x = Id(self,g,x)
            if g.supp(1) <= x && x<= g.supp(2)
                Ig_x = g.eval(x)/2;
            else
                Ig_x = 0;
            end
        end
        
        function f = get_RHS_data(self,uinc)
            f = CombinedFunction(uinc,self.domain,self.coupling_param);
        end
        
        function X = getSymmetries(self)
            X = self.domain.getSymmetries();
            d = abs(min(min(X)))+max(max(X))+1;
            X(X<0) = X(X<0)+d;
        end
    end
end