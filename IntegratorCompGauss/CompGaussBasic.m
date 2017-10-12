
classdef CompGaussBasic < Integrator
    %Composite Gauss integrator, basic and inefficient, for testing
    
    properties
        wavelength   %will want to increase points per wavelength
        pointsPerWavelength
        pointsPerWavelengthSmooth
        memLim=2000000
    end
    
    methods
        function self=CompGaussBasic(kwave,pointsPerWavelength,pointsPerWavelengthSmooth)
            if nargin==2
                self.pointsPerWavelengthSmooth=pointsPerWavelength;
            else
                self.pointsPerWavelengthSmooth=pointsPerWavelengthSmooth;
            end
            self.pointsPerWavelength=pointsPerWavelength;
            self.wavelength=2*pi/kwave;            
        end
        
        function val=eval(self,integral, X)
            if isa(integral,'BoundaryIntegral')
                %ideally would look at properties of integral here...

                %check for if parametrised values have been given:
                if integral.toR2==0
                   X=integral.domain{end}.trace(X); 
                end

                supp=integral.boundaryFn.supp;

                %now create weights and nodes:
                [s, w]=gauss_quad_wave_split(supp(1), supp(2), self.wavelength, self.pointsPerWavelength );

                 Y_s=integral.domain.trace(s);

                [X_1,Y_1]=meshgrid(X(:,1), Y_s(:,1));
                [X_2,Y_2]=meshgrid(X(:,2), Y_s(:,2));

                R=sqrt( (X_1-Y_1).^2 + (X_2-Y_2).^2 );
                val=integral.kernel(R')*(w.*integral.boundaryFn.eval(s));
            elseif isa(integral,'BEMintegral3D') && isequal(integral.domain{1},integral.Fn2.domain)
                supp_s=integral.Fn1.supp;
                supp_q=integral.domain{2}.supp;
                supp_t=integral.Fn2.supp;
                if integral.singularKerIndex==1
                    Ns=self.pointsPerWavelength;
                    Nt=self.pointsPerWavelengthSmooth;                    
                elseif integral.singularKerIndex==2
                    Nt=self.pointsPerWavelength;
                    Ns=self.pointsPerWavelengthSmooth;   
                end
                [s, w_s]=gauss_quad_wave_split(supp_s(1), supp_s(2), self.wavelength, Ns, integral.Fn1.suppWidth );
                [q, w_q]=gauss_quad_wave_split(supp_q(1), supp_q(2), self.wavelength, self.pointsPerWavelength +1, integral.Fn1.suppWidth );
                [t, w_t]=gauss_quad_wave_split(supp_t(1), supp_t(2), self.wavelength, Nt, integral.Fn2.suppWidth);
                [S, Q, T] = kron3( s,q,t ); [wS, wQ, wT] = kron3( w_s, w_q, w_t );
                if length(S)>self.memLim
                    val=0;
                    splits=ceil(length(S)/self.memLim);
                    splitLength=ceil(length(S)/splits);
                    for j=1:(splits-1)
                        indices=(1+(j-1)*splitLength):(j*splitLength);
                        val=val+(wS(indices).*wQ(indices).*wT(indices)).'*(integral.kernel{1}.evalFull(S(indices),Q(indices),integral.domain{1},integral.domain{2}).*integral.kernel{2}.evalFull(Q(indices),T(indices),integral.domain{2},integral.domain{3}).*integral.Fn1.eval(S(indices)).*conj(integral.Fn2.eval(T(indices))));
                    end
                    indices=(1+(splits-1)*splitLength):length(S);
                    val=val+(wS(indices).*wQ(indices).*wT(indices)).'*(integral.kernel{1}.evalFull(S(indices),Q(indices),integral.domain{1},integral.domain{2}).*integral.kernel{2}.evalFull(Q(indices),T(indices),integral.domain{2},integral.domain{3}).*integral.Fn1.eval(S(indices)).*conj(integral.Fn2.eval(T(indices))));
                else
                    val=(wS.*wQ.*wT).'*(integral.kernel{1}.evalFull(S,Q,integral.domain{1},integral.domain{2}).*integral.kernel{2}.evalFull(Q,T,integral.domain{2},integral.domain{3}).*integral.Fn1.eval(S).*conj(integral.Fn2.eval(T)));
                end
            elseif isa(integral,'BEMintegral2D') && isequal(integral.domain{1},integral.Fn2.domain)
                supp_s=integral.Fn1.supp;
                supp_t=integral.Fn2.supp;
                [s, w_s]=gauss_quad_wave_split(supp_s(1), supp_s(2), self.wavelength, self.pointsPerWavelength, integral.Fn1.suppWidth );
                [t, w_t]=gauss_quad_wave_split(supp_t(1), supp_t(2), self.wavelength, self.pointsPerWavelength +1, integral.Fn2.suppWidth);
                [S,T] = kron2( s,t ); [wS,wT] = kron2( w_s,w_t );
                val=(wS.*wT).'*(integral.kernel.evalFull(S,T,integral.domain{1},integral.domain{2}).*integral.Fn1.eval(S).*conj(integral.Fn2.eval(T)));
            elseif isa(integral,'InnerProduct1D')
                if length(integral.supp)==2 %check that integration domain has positive measure
                    [s, w]=gauss_quad_wave_split(integral.supp(1), integral.supp(2), self.wavelength, self.pointsPerWavelength );
                    val=w.'*(integral.Fn1.eval(s).*conj(integral.Fn2.eval(s)));
                else
                    val=0;
                end
            elseif isa(integral,'BEMintegral1D')
                if isequal(integral.SubIntegrals,[])
                    [ x, w ]=gauss_quad_wave_split(integral.supp(1), integral.supp(2), self.wavelength, self.pointsPerWavelength, integral.suppWidth );
                else
                    [ x1, w1 ]=gauss_quad_wave_split(integral.SubIntegrals{1}.supp(1), integral.SubIntegrals{1}.supp(2), self.wavelength, self.pointsPerWavelength, integral.SubIntegrals{1}.suppWidth );
                    [ x2, w2 ]=gauss_quad_wave_split(integral.SubIntegrals{2}.supp(1), integral.SubIntegrals{2}.supp(2), self.wavelength, self.pointsPerWavelength, integral.SubIntegrals{2}.suppWidth );
                    x=[x1; x2;]; w=[w1; w2;];
                end
                
                val=w.'*(integral.kernel.evalFull(integral.ColPoint,x,integral.domain,integral.domain).*integral.Fn.eval(x));
            else
                error('No routine seems to exist for integral');
            end
        end
    end
    
end