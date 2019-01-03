classdef  baseFnHNA < BoundaryFunction
    properties
        p
        kwave
        %suppWidth      %in parent class
        pm      %essentially the phase function, can be {1,-1,0}
        normaliser
        phaseLinear
        L       %original width of Legendre basis element
    end
    
    methods 
        function self=baseFnHNA(kwave,p,meshEl,pm,side)
            self.p=p;
            self.kwave=kwave;
            self.pm=pm;
            self.phaseLinear=[pm 0];
            self.supp=meshEl.interval;
            self.suppWidth=meshEl.width;
            %now store the supp data as seperate varaibales, incase this
            %function is later restricted to a subdomain; we don't want the
            %value of the function (given my 'eval' below) to change
            self.a=self.supp(1); self.b=self.supp(2); self.L=self.suppWidth;
            self.normaliser=sqrt((2*p + 1)/2);
            self.oscillator=@(s) exp(1i*kwave*pm*s);
            self.domain=side;
            self.phaseMaxStationaryPointOrder=0;
            self.meshEl = meshEl;
            %self.phase = {@(x) pm*x, @(x) pm, @(x) 0};
        end
        
        function [y, yNonOsc]=eval(self,x)
           yNonOsc=self.normaliser*sqrt(2/self.L)*(x>=self.supp(1)).*(x<=self.supp(2)).*legendre(self.p, 2*((x-self.a)./(self.L)) -1);
           y=yNonOsc.*exp(self.pm*1i*self.kwave.*x) ;
        end
        
        function y = evalNonOscAnalPivot(self,z,~,pivLminus_a)
            %if b-a is small, but 0<<a, then rounding errors can occur.
            %This function essentially returns
            %evalNonOscAnal(self,pivL-x,~), in the case where pivL-x rounds
            %to zero, avoiding the rounding error
            newArg = pivLminus_a/self.L - z/self.L;
            y = self.normaliser*sqrt(2/self.L)*legendre(self.p, 2*(newArg) -1);
        end
        
        function y = evalNonOscAnal(self,x,~)
            y = self.normaliser*sqrt(2/self.L)*legendre(self.p, 2*((x-self.a)./(self.L)) -1);
        end
        
        function y = evalNonOscAnalPivot(self,z,~,pivLminus_a)
            %if b-a is small, but 0<<a, then rounding errors can occur.
            %This function essentially returns
            %evalNonOscAnal(self,pivL-x,~), in the case where pivL-x rounds
            %to zero, avoiding the rounding error
            newArg = pivLminus_a/self.L - z/self.L;
            y = self.normaliser*sqrt(2/self.L)*legendre(self.p, 2*(newArg) -1);
        end
        
        function g = phaseAnal(self,s,deriv,~)
            switch deriv
                case 0
                    g = self.pm*s;
                case 1
                    g = self.pm;
                otherwise
                    g = 0;
            end
        end
    end
    
end

