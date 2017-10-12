classdef BEMintegral2D < BEMintegral
    %decomposes integral into nice(r) components, to be fed into an
    %'integrator' type object
    
    properties 
        subIntegrals
        nearParam=0.15
    end
    
    methods
        function self=BEMintegral2D(intOp, basFun1, basFun2, integralType)
            %classify what type of integral this is, and split into
            %subintegrals if we need, should be at most one recursion
                
            if nargin<=3
                self.type='unknown';
            else
                self.type=integralType;
            end
            
            self.domain{1}=basFun1.supp;
            self.domain{2}=basFun2.supp;
            self.kernel=intOp.kernel;

            if strcmp(self.type,'unknown');
                %if first and last sides, make an extra side to avoid all the
                %extra weird cases
                %if basFun1.side-basFun2.side
%                    warning('remember to add extra side for first and last side for polygons');
                %end

                splitPoints = domainSplitter( self.domain );

                splits1=length(splitPoints{1})-1;
                splits2=length(splitPoints{2})-1;

                if splits1==1 && splits2==1
                    self.type=integralClassify(self.domain{1},self.domain{2}, self.nearParam);
                else
                    %loop over the subdomains, make the subIntegrals
                    subIntCounter=0;
                    for n=1:splits1
                        for m=1:splits2
                            subIntCounter=subIntCounter+1;
                            subDom1=[splitPoints{1}(n) splitPoints{1}(n+1)];
                            subDom2=[splitPoints{2}(m) splitPoints{2}(m+1)];
                            subType=integralClassify(subDom1,subDom2, self.nearParam);
                            self.subIntegrals{subIntCounter}=BEMintegral2D(intOp, restrictTo(basFun1,subDom1), restrictTo(basFun2,subDom2), subType);
                        end
                    end
                end
            end
            
            
        end
    end
    
end

