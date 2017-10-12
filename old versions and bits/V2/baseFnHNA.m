classdef  baseFnHNA < BoundaryFunction
    properties
        p
        kwave
        %suppWidth      %in parent class
        pm      %essentially the phase function, can be {1,-1,0}
        normaliser
    end
    
    methods 
        function self=baseFnHNA(kwave,p,meshEl,pm,side)
            self.p=p;
            self.kwave=kwave;
            self.pm=pm;
            self.supp=meshEl.interval;
            self.suppWidth=meshEl.width;
            self.normaliser=sqrt((2*p + 1)/2);
            self.oscillator=@(s) exp(1i*kwave*pm*s);
            self.domain=side;
        end
        
        function y=eval(self,x)
           %val=legendre(p,s).*exp(pm*1i*self.kwave.*s) ;
           y=self.normaliser*sqrt(2/self.suppWidth)*(x>=self.supp(1)).*(x<=self.supp(2)).*legendre(self.p, 2*((x-self.supp(1))./(self.suppWidth)) -1).*exp(self.pm*1i*self.kwave.*x) ;
        end
    end
    
end

