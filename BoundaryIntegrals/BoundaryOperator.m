classdef BoundaryOperator < handle
    %general BIE class. Not quite abstract.
    %NEED TO REWRITE SOME OF THE COMMENTED BITS
    
    properties
        kernel
        domain
    end
    
    methods
        function self=BoundaryOperator(kernel_,domain)
             if length(domain)==1
                 self.domain{1}=domain;
                 self.domain{2}=domain;
             else 
                 self.domain=domain;
             end
             self.kernel=kernel_;
             %self.kwave=kwave;
        end
        
        function K=mtimes(X,Y)
            %need a check that 'c' is 1D here
            if isequal(size(X),[1 1]) && isa(X,'numeric')
                K=BoundaryOperator(X*Y.kernel,Y.domain);
%                 kernelScaled=@(kwave,x,y,extras) c*A.kernel(kwave,x,y,extras);
%                 K=BoundaryOperator(A.kwave,kernelScaled);
            elseif isa(X,'BoundaryOperator') && isa(Y,'BoundaryOperator')
                if isequal (X.domain{end},Y.domain{1}) 
                    K=BoundaryIntegral(X.kernel,Y,X.domain);
                else
                    error('Domain of BoundaryFunction must match operator');
                end
                %compatible operators, but need to check 'isequal' works OK
                %here...
                timesDomain=X;
                for n=(length(X)+1):(length(Y)+length(X)-1)
                    timesDomain{n}=Y(n-length(X)+1);
                end
                timesKernel=X.kernel;
                for n=1:length(Y.kernel)
                    timesKernel{n+length(Y.kernel)}=Y.kernel{n};
                end
                %now put it all together
                K=BoundaryOperator(timesKernel,timesDomain);
            elseif isa(X,'BoundaryOperator') && isa(Y,'BoundaryFunction')
                if isequal(Y.domain,X.domain{end})
                    K=BoundaryIntegral(X.kernel,Y,X.domain);
                else
                    error('Domain of BoundaryFunction must match operator');
                end
            elseif isa(X,'BoundaryOperator') && isa(Y,'BoundaryIntegral')
                %error('Havent coded this yet');
                if isequal(X.domain{end},Y.domain{1})
                    if isa(X.kernel,'cell')
                        for nx=1:length(X.kernel)
                            newKernel{nx}=X.kernel{nx};                        
                        end
                    else
                        newKernel{1}=X.kernel;
                    end
                    if isa(Y.kernel,'cell')
                        nyCount=1;
                        for ny=(length(X.kernel)+1):(length(X.kernel)+length(Y.kernel))
                            newKernel{ny}=Y.kernel{nyCount};
                            nyCount=nyCount+1;
                        end
                    else
                        newKernel{length(newKernel)+1}=Y.kernel;
                    end
                    
                    if isa(X.domain,'cell')
                        for nx=1:(length(X.domain)-1)
                            newDomain{nx}=X.domain{nx};                        
                        end
                    else
                        newDomain{1}=X.domain;                        
                    end
                    if isa(Y.domain,'cell')
                        nyCount=1;
                        for ny=(length(X.domain)):(length(X.domain)+length(Y.domain)-1)
                            newDomain{ny}=Y.domain{nyCount};           
                            nyCount=nyCount+1;
                        end           
                    else
                        newDomain{length(newDomain)+1}=Y.kernel;
                    end
                    K=BoundaryIntegral(newKernel,Y.Fn,newDomain);
                else
                    error('Domain of BoundaryFunction must match operator');
                end
            else
                error('cannot multiply operator by vector');
            end
            
        end
    end
    
end

