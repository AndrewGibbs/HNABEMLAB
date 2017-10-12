classdef kernel
    %kernel class, designed to be used in conjunction with BoundaryOperator
    
    properties %all properties will be handle functions
        S    %same side
        D    %disjoint side
        N    %neighbouring side
        R    %radial R2 version
    end
    
    methods
        function self=kernel(S,D,N,R)
            self.S=S;
            self.D=D;
            self.N=N;
            self.R=R;
        end
        
        function kernelSum=plus(self1,self2)
            %S means same side kernel
            sumS=@(s,t,r) self1.kerS(s,t,r)+self2.kerS(s,t,r);
            %D means disjoint sides kernel
            sumD=@(s,t,Gs,Gt) self1.kerS(s,t,Gs,Gt)+self2.kerS(s,t,Gs,Gt);
            %N means neighbouring sides kernel
            sumN=@(s,t,ns,dSdt) self1.kerS(s,t,ns,dSdt)+self2.kerS(s,t,ns,dSdt);
            sumR=@(r) self1.kerR(s,t,r)+self2.kerR(r);
            %create new kernel
            kernelSum=kernel(sumS,sumD,sumN,sumR);
        end
        
        
        function kernelTimes=mtimes(c,self)
            if isequa(size(c),[1,1])
                cS=@(s,t,r) c*self.kerS(s,t,r);
                cD=@(s,t,Gs,Gt) c*self.kerS(s,t,Gs,Gt);
                cN=@(s,t,ns,dSdt) c*self.kerS(s,t,ns,dSdt);
                cR=@(r) c*self.kerR(r);
                kernelTimes=kernel(cS,cD,cN,cR);
            else
                error('scalar multiplication of kernels only');
            end
        end
    end
    
end