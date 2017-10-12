function out=SingleLayer(kwave,domain)
    %single layer operator defined on a single domain

    %self.kernel=@(x,y) (1i/4)*besselh(0,1,kwave*abs(x-y));
    ker=SLkernel(kwave);
    out= BoundaryOperator( ker,domain);
    %self.domain=domain;

    
    
end

