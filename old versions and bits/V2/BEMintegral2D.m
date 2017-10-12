classdef BEMintegral2D < BEMintegral
    %decomposes integral into nice(r) components, to be fed into an
    %'integrator' type object
    
    properties 
        subIntegrals
        nearParam=0.15
        Fn1
        Fn2
        supp
        suppWidth
    end
    
    methods
        function self=BEMintegral2D(domain, kernel, Fn1, Fn2, integralType)
            %classify what type of integral this is, and split into
            %subintegrals if we need, should be at most one recursion
                
            if nargin<=4
                self.type='unknown';
            else
                self.type=integralType;
            end
            
            self.supp{1}=Fn1.supp;
            self.supp{2}=Fn2.supp;
            self.domain=domain;
            if (~isequal(Fn1.domain,domain{1})) || (~isequal(Fn2.domain,domain{2}))
                error('Domains do not match functions');
            end
            self.Fn1=Fn1;
            self.Fn2=Fn2;
            self.suppWidth=[Fn1.suppWidth Fn2.suppWidth];
            self.kernel=kernel;
            %need to extract kernel from oscillator too!
            self.oscillator=@(s,t) self.Fn1.oscillator(s)-self.Fn2.oscillator(t) ;

            if strcmp(self.type,'unknown');
                %if first and last sides, make an extra side to avoid all the
                %extra weird cases
                %if basFun1.side-basFun2.side
%                    warning('remember to add extra side for first and last side for polygons');
                %end

                splitPoints = domainSplitter( self.supp );

                splits1=length(splitPoints{1})-1;
                splits2=length(splitPoints{2})-1;

                if splits1==1 && splits2==1
                    self.type=integralClassify(self.supp{1},self.supp{2}, self.nearParam);
                else
                    %loop over the subdomains, make the subIntegrals
                    subIntCounter=0;
                    for n=1:splits1
                        for m=1:splits2
                            subIntCounter=subIntCounter+1;
                            subDom1=[splitPoints{1}(n) splitPoints{1}(n+1)];
                            subDom2=[splitPoints{2}(m) splitPoints{2}(m+1)];
                            subType=integralClassify(subDom1,subDom2, self.nearParam);
                            self.subIntegrals{subIntCounter}=BEMintegral2D(domain, kernel, restrictTo(Fn1,subDom1), restrictTo(Fn2,subDom2), subType);
                        end
                    end
                end
            end
            
            
        end
        
    end
    
end

