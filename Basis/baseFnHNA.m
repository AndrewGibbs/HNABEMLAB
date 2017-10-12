classdef  baseFnHNA < BoundaryFunction
    properties
        p
        kwave
        %suppWidth      %in parent class
        pm      %essentially the phase function, can be {1,-1,0}
        normaliser
        phaseLinear
        a       %left-hand endpoints of the support of the Legendre bsais element
        b       %right-hand endpoints of the support of the Legendre bsais element
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
        end
        
        function [y, yNonOsc]=eval(self,x)
           %val=legendre(p,s).*exp(pm*1i*self.kwave.*s) ;
           %yNonOsc=self.normaliser*sqrt(2/self.L)*(x>=self.supp(1)).*(x<=self.supp(2)).*legendre(self.p, 2*((x-self.a)./(self.L)) -1);
           yNonOsc=self.normaliser*sqrt(2/self.L)*(x>=self.supp(1)).*(x<=self.supp(2)).*legendre(self.p, 2*((x-self.a)./(self.L)) -1);
           y=yNonOsc.*exp(self.pm*1i*self.kwave.*x) ;
           %y=self.normaliser*sqrt(2/self.L)*(x>=self.supp(1)).*(x<=self.supp(2)).*legendre(self.p, 2*((x-self.a)./(self.L)) -1).*exp(self.pm*1i*self.kwave.*x) ;
        end
    end
    
end

